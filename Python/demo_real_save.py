import argparse
import os
import time
import traceback
from pathlib import Path
import numpy as np
import imageio.v2 as imageio
from scipy.io import savemat

# Try to import CuPy.
try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False

from utils.load_cpfa import load_cpfa
from utils.init_interp import init_interp
from para_set import para_set
from utils.init_dict import init_dict
from utils.image2patch import image2patch
from utils.NLCSR import NLCSR
from utils.patch2image import patch2image
from utils.S_update import S_update
from utils.cal_Stokes import cal_Stokes
from utils.matlab_compat import im2uint8
from utils.colorjetmap import colorjetmap
from utils.contrast_str import contrast_str


def _imwrite(img, path):
    if img.dtype != np.uint8:
        img = im2uint8(img)
    # Ensure that the destination parent directory exists.
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    imageio.imwrite(path, img)


def _to_numpy(obj):
    if _HAS_CUPY and isinstance(obj, cp.ndarray):
        return cp.asnumpy(obj)
    if isinstance(obj, dict):
        return {k: _to_numpy(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        vals = [_to_numpy(v) for v in obj]
        return type(obj)(vals)
    return obj


def main():
    parser = argparse.ArgumentParser(description='real demo (DFT-only, optional CUDA)')
    parser.add_argument('--indices', type=str, default='0', help='comma-separated image indices')
    parser.add_argument('--image-names', type=str, default='', help='comma-separated image names')
    parser.add_argument('--image-root', type=str, default=os.path.join('data', 'real'), help='root folder')
    parser.add_argument('--gpu-id', type=int, default=None, help='worker GPU id')
    parser.add_argument('--use-cuda', action='store_true', help='use CuPy')
    parser.add_argument('--save-dir', type=str, default='real_resultsPython_cuda_dft')
    args = parser.parse_args()

    os.environ['USE_EXPLICIT_DFT'] = '1'
    os.environ['USE_MATLAB_ENGINE_IFFT'] = '0'
    os.environ.setdefault('CUPY_CACHE_DIR', str((Path(__file__).resolve().parent / '.cupy_cache')))

    if args.gpu_id is not None:
        os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu_id)
    
    if args.use_cuda:
        if not _HAS_CUPY:
            print("[Warning] --use-cuda requested but cupy not found. Falling back to CPU.")
            os.environ['USE_CUDA'] = '0'
            use_cuda = False
        else:
            os.environ['USE_CUDA'] = '1'
            use_cuda = True
    else:
        os.environ['USE_CUDA'] = '0'
        use_cuda = False

    # save_dir is the output root (or a temporary root assigned by the scheduler).
    save_dir_root = Path(args.save_dir)
    save_dir_root.mkdir(parents=True, exist_ok=True)

    image_root = Path(args.image_root)
    if not image_root.exists():
        raise ValueError(f'image root not found: {image_root}')

    image_paths = sorted([p for p in image_root.rglob('*.png') if p.is_file()])
    if not image_paths:
        raise ValueError(f'no png images found under: {image_root}')

    image_rel = [str(p.relative_to(image_root)).replace('\\', '/') for p in image_paths]

    if args.image_names.strip():
        req = [x.strip() for x in args.image_names.split(',') if x.strip() != '']
        norm_rel = {n.lower(): i for i, n in enumerate(image_rel)}
        norm_base = {}
        for i, n in enumerate(image_rel):
            b = Path(n).name.lower()
            norm_base.setdefault(b, []).append(i)
        indices = []
        for name in req:
            key = name.lower()
            if not key.endswith('.png'):
                key = key + '.png'
            if key in norm_rel:
                indices.append(norm_rel[key])
                continue
            cand = norm_base.get(Path(key).name, [])
            if len(cand) == 1:
                indices.append(cand[0])
                continue
            if len(cand) > 1:
                raise ValueError(f'ambiguous image name "{name}", please use relative path under image root')
            raise ValueError(f'image not found under {image_root}: {name}')
    else:
        indices = [int(x) for x in args.indices.split(',') if x.strip() != '']

    t0 = time.perf_counter()

    for nn in indices:
        try:
            t_img = time.perf_counter()
            image_name = image_rel[nn]
            
            # 1. Load and initialize.
            cpfa = load_cpfa(str(image_root / image_name))
            s_ini, mask = init_interp(cpfa)
            para = para_set(s_ini)
            
            if use_cuda:
                cpfa_g = cp.asarray(cpfa, dtype=cp.float64)
                mask_g = cp.asarray(mask, dtype=cp.float64)
                s_ini_g = cp.asarray(s_ini, dtype=cp.float64)
                para.Smask = cpfa_g[:, :, None] * mask_g
                mask_used = mask_g
                s = s_ini_g
            else:
                para.Smask = cpfa[:, :, None] * mask
                mask_used = mask
                s = s_ini
                
            d0 = init_dict(para)
            d = cp.asarray(d0, dtype=cp.float64) if use_cuda else d0

            # 2. Iterative reconstruction.
            for _ in range(para.main_iternum):
                s_pat, para = image2patch(s, para)
                if use_cuda and not isinstance(s_pat, cp.ndarray):
                    s_pat = cp.asarray(s_pat)

                l = s_pat.shape[3]
                patch_num = para.patch_num if l >= 1000 else l
                n_patches = int(np.ceil(l / patch_num))
                
                if use_cuda:
                    d = cp.tile(d[..., cp.newaxis], (1, 1, 1, n_patches))
                    dx = cp.zeros_like(s_pat)
                else:
                    d = np.tile(d[..., np.newaxis], (1, 1, 1, n_patches))
                    dx = np.zeros_like(s_pat)
                
                for n in range(n_patches):
                    idx_start = n * patch_num
                    idx_end = min((n + 1) * patch_num, l)
                    dx_cur, d_cur = NLCSR(
                        d[:, :, :, n],
                        s_pat[:, :, :, idx_start:idx_end],
                        para,
                        return_numpy=(not use_cuda),
                    )
                    dx[:, :, :, idx_start:idx_end] = dx_cur
                    d[:, :, :, n] = d_cur
                    
                s_dx = patch2image(dx, para)
                s_rec = S_update(s_dx, mask_used, para)
                s = s_rec

            # 3. Transfer results back to the CPU.
            if use_cuda:
                cp.cuda.Stream.null.synchronize()
                s_rec_cpu = cp.asnumpy(s_rec)
                s_ini_cpu = cp.asnumpy(s_ini_g)
            else:
                s_rec_cpu = s_rec
                s_ini_cpu = s_ini

            s0, dolp, aolp = cal_Stokes(s_rec_cpu)
            ini_s0, ini_dolp, ini_aolp = cal_Stokes(s_ini_cpu)

            # ==========================================
            # Preserve the input subdirectory structure in the output path.
            # ==========================================
            p_rel = Path(image_name) # e.g., CPDnet_input001/1.png
            subdir = p_rel.parent    # e.g., CPDnet_input001
            stem = p_rel.stem        # e.g., 1
            
            # Create the corresponding subdirectory under save_dir.
            current_out_dir = save_dir_root / subdir
            current_out_dir.mkdir(parents=True, exist_ok=True)
            
            # Use only the image name as the base filename (for example, '1').
            base = stem 
            
            # Construct the complete output path.
            _imwrite(s_rec_cpu[:, :, 6:9], str(current_out_dir / f'{base}_90.png'))
            _imwrite(s_rec_cpu[:, :, 3:6], str(current_out_dir / f'{base}_45.png'))
            _imwrite(s_rec_cpu[:, :, 9:12], str(current_out_dir / f'{base}_135.png'))
            _imwrite(s_rec_cpu[:, :, 0:3], str(current_out_dir / f'{base}_0.png'))
            _imwrite(s0, str(current_out_dir / f'{base}_S0.png'))
            _imwrite(colorjetmap(contrast_str(dolp, 0.3)), str(current_out_dir / f'{base}_DoLP.png'))
            _imwrite(colorjetmap((aolp + np.pi / 2) / np.pi), str(current_out_dir / f'{base}_AoLP.png'))

            savemat(
                str(current_out_dir / f'{base}_intermediate.mat'),
                {
                    'CPFA': cpfa,
                    'Mask': mask,
                    'S_ini': s_ini_cpu,
                    'S_rec': s_rec_cpu,
                    'S0': s0,
                    'DoLP': dolp,
                    'AoLP': aolp,
                    'iniS0': ini_s0,
                    'iniDoLP': ini_dolp,
                    'iniAoLP': ini_aolp,
                    'para': _to_numpy(para.__dict__),
                },
            )
            print(f'image={image_name}, elapsed_sec={time.perf_counter() - t_img:.3f}')

        except KeyboardInterrupt:
            raise
        except Exception as e:
            print(f"[ERROR] Failed to process image: {image_name}")
            print(f"        Reason: {str(e)}")
            traceback.print_exc()
            continue

    elapsed = time.perf_counter() - t0
    print(f'total_elapsed_sec={elapsed:.3f}')


if __name__ == '__main__':
    main()
