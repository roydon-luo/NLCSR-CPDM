import argparse
import os
import time
import traceback
from pathlib import Path

import imageio.v2 as imageio
import numpy as np
# from scipy.io import savemat  # Intentionally disabled.

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False

from para_set import para_set
from utils.NLCSR import NLCSR
from utils.S_update import S_update
from utils.colorjetmap import colorjetmap
from utils.contrast_str import contrast_str
from utils.image2patch import image2patch
from utils.init_dict import init_dict
from utils.matlab_compat import im2double, im2uint8, rgb2gray
from utils.patch2image import patch2image


def _imwrite(img, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if img.dtype != np.uint8:
        img = im2uint8(img)
    imageio.imwrite(path, img)


def _mat2gray(img):
    img = np.asarray(img, dtype=np.float64)
    vmin = float(np.min(img))
    vmax = float(np.max(img))
    if vmax <= vmin:
        return np.zeros_like(img, dtype=np.float64)
    return (img - vmin) / (vmax - vmin)


def _collect_images(image_root: Path):
    exts = {".bmp", ".png", ".jpg", ".jpeg", ".tif", ".tiff"}
    paths = [p for p in image_root.rglob("*") if p.is_file() and p.suffix.lower() in exts]
    return sorted(paths)


def _calc_stokes_gray(s4):
    # channel order: [I0, I45, I90, I135]
    i0 = s4[:, :, 0]
    i45 = s4[:, :, 1]
    i90 = s4[:, :, 2]
    i135 = s4[:, :, 3]
    s0 = 0.5 * (i0 + i45 + i90 + i135)
    s1 = i0 - i90
    s2 = i45 - i135
    dop = np.sqrt(s1**2 + s2**2) / (s0 + np.finfo(np.float64).eps)
    aop = np.mod(0.5 * np.arctan2(s2, s1), np.pi)
    return s0, dop, aop


def main():
    parser = argparse.ArgumentParser(description="real gray polarization denoise (no interpolation)")
    parser.add_argument(
        "--image-root",
        type=str,
        default=os.path.join("data", "Noise", "temp-09172025123715-0000"),
        help="root folder of raw polarization images",
    )
    # Note: this is the singular --gpu-id option, not --gpus.
    parser.add_argument("--indices", type=str, default="", help="comma-separated indices (optional)")
    parser.add_argument("--image-names", type=str, default="", help="comma-separated image names (optional)")
    parser.add_argument("--gpu-id", type=int, default=None, help="worker GPU id")
    parser.add_argument("--use-cuda", action="store_true", help="use CuPy")
    parser.add_argument("--save-dir", type=str, default="real_resultsPython_grayDenoise_nointerp")
    args = parser.parse_args()

    os.environ["USE_EXPLICIT_DFT"] = "1"
    os.environ["USE_MATLAB_ENGINE_IFFT"] = "0"
    os.environ.setdefault("CUPY_CACHE_DIR", str((Path(__file__).resolve().parent / ".cupy_cache")))

    if args.gpu_id is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

    if args.use_cuda:
        if not _HAS_CUPY:
            print("[Warning] --use-cuda requested but cupy not found. Falling back to CPU.")
            os.environ["USE_CUDA"] = "0"
            use_cuda = False
        else:
            os.environ["USE_CUDA"] = "1"
            use_cuda = True
    else:
        os.environ["USE_CUDA"] = "0"
        use_cuda = False

    image_root = Path(args.image_root)
    if not image_root.exists():
        raise ValueError(f"image root not found: {image_root}")

    save_root = Path(args.save_dir)
    save_root.mkdir(parents=True, exist_ok=True)

    image_paths = _collect_images(image_root)
    if not image_paths:
        raise ValueError(f"no image files found under: {image_root}")

    image_rel = [str(p.relative_to(image_root)).replace("\\", "/") for p in image_paths]

    if args.image_names.strip():
        req = [x.strip() for x in args.image_names.split(",") if x.strip()]
        norm_rel = {n.lower(): i for i, n in enumerate(image_rel)}
        norm_base = {}
        for i, n in enumerate(image_rel):
            b = Path(n).name.lower()
            norm_base.setdefault(b, []).append(i)
        indices = []
        for name in req:
            key = name.lower()
            if key in norm_rel:
                indices.append(norm_rel[key])
                continue
            cand = norm_base.get(Path(key).name, [])
            if len(cand) == 1:
                indices.append(cand[0])
                continue
            if len(cand) > 1:
                raise ValueError(f'ambiguous image name "{name}", please use relative path under image root')
            raise ValueError(f"image not found under {image_root}: {name}")
    elif args.indices.strip():
        indices = [int(x) for x in args.indices.split(",") if x.strip()]
    else:
        indices = list(range(len(image_paths)))

    t0 = time.perf_counter()
    skipped_count = 0
    processed_count = 0

    for nn in indices:
        image_name = image_rel[nn]
        t_img = time.perf_counter()
        try:
            p_rel = Path(image_name)
            subdir = p_rel.parent
            stem = p_rel.stem
            suffix = p_rel.suffix if p_rel.suffix else ".png"

            dir_raw = save_root / "Processed_Raw" / subdir
            dir_i0 = save_root / "Split_I0" / subdir
            dir_i45 = save_root / "Split_I45" / subdir
            dir_i90 = save_root / "Split_I90" / subdir
            dir_i135 = save_root / "Split_I135" / subdir
            dir_aop = save_root / "Phys_AoP" / subdir
            dir_dop = save_root / "Phys_DoP" / subdir
            
            # Ensure that all output directories exist.
            for d in [dir_raw, dir_i0, dir_i45, dir_i90, dir_i135, dir_aop, dir_dop]:
                d.mkdir(parents=True, exist_ok=True)

            target_check = dir_raw / f"{stem}{suffix}"
            if target_check.exists():
                skipped_count += 1
                print(f"[Skip] image={image_name}, reason=already_exists")
                continue

            raw = imageio.imread(str(image_root / image_name))
            raw = im2double(raw)
            if raw.ndim == 3:
                raw = rgb2gray(raw)

            h, w = raw.shape
            if (h % 2) != 0 or (w % 2) != 0:
                print(f"[Skip] image={image_name}, reason=odd_size({h}x{w})")
                skipped_count += 1
                continue

            # Sony IMX250MZR pattern:
            # [90, 45]
            # [135, 0]
            i90 = raw[0::2, 0::2]
            i45 = raw[0::2, 1::2]
            i135 = raw[1::2, 0::2]
            i0 = raw[1::2, 1::2]
            s_split = np.stack([i0, i45, i90, i135], axis=2)

            para = para_set(s_split)
            para.Smask = s_split
            para.rho_hub = np.repeat(0.05, para.ch)
            para.delta_hub = np.repeat(0.001, para.ch)
            para.beta_hub = np.repeat(0.95, para.ch)

            if use_cuda:
                s = cp.asarray(s_split, dtype=cp.float64)
                mask_used = cp.ones_like(s)
                para.Smask = cp.asarray(s_split, dtype=cp.float64)
                d = cp.asarray(init_dict(para), dtype=cp.float64)
            else:
                s = s_split
                mask_used = np.ones_like(s_split)
                d = init_dict(para)

            for _ in range(para.main_iternum):
                s_pat, para = image2patch(s, para)
                if use_cuda and not isinstance(s_pat, cp.ndarray):
                    s_pat = cp.asarray(s_pat)
                l = s_pat.shape[3]
                patch_num = para.patch_num if l >= 1000 else l
                n_patches = int(np.ceil(l / patch_num))

                if use_cuda:
                    d_ext = cp.tile(d[..., cp.newaxis], (1, 1, 1, n_patches))
                    dx = cp.zeros_like(s_pat)
                else:
                    d_ext = np.tile(d[..., np.newaxis], (1, 1, 1, n_patches))
                    dx = np.zeros_like(s_pat)

                for n in range(n_patches):
                    idx_start = n * patch_num
                    idx_end = min((n + 1) * patch_num, l)
                    dx_cur, d_cur = NLCSR(
                        d_ext[:, :, :, n],
                        s_pat[:, :, :, idx_start:idx_end],
                        para,
                        return_numpy=(not use_cuda),
                    )
                    dx[:, :, :, idx_start:idx_end] = dx_cur
                    d_ext[:, :, :, n] = d_cur

                s_dx = patch2image(dx, para)
                s_rec_sub = S_update(s_dx, mask_used, para)
                s = s_rec_sub

            if use_cuda:
                cp.cuda.Stream.null.synchronize()
                s_rec_cpu = cp.asnumpy(s_rec_sub)
            else:
                s_rec_cpu = s_rec_sub

            denoised_raw = np.zeros((h, w), dtype=np.float64)
            denoised_raw[0::2, 0::2] = s_rec_cpu[:, :, 2]  # 90
            denoised_raw[0::2, 1::2] = s_rec_cpu[:, :, 1]  # 45
            denoised_raw[1::2, 0::2] = s_rec_cpu[:, :, 3]  # 135
            denoised_raw[1::2, 1::2] = s_rec_cpu[:, :, 0]  # 0

            s0, dop, aop = _calc_stokes_gray(s_rec_cpu)

            _imwrite(_mat2gray(denoised_raw), str(target_check))
            _imwrite(_mat2gray(s_rec_cpu[:, :, 0]), str(dir_i0 / f"{stem}.png"))
            _imwrite(_mat2gray(s_rec_cpu[:, :, 1]), str(dir_i45 / f"{stem}.png"))
            _imwrite(_mat2gray(s_rec_cpu[:, :, 2]), str(dir_i90 / f"{stem}.png"))
            _imwrite(_mat2gray(s_rec_cpu[:, :, 3]), str(dir_i135 / f"{stem}.png"))
            _imwrite(colorjetmap(contrast_str(dop, 1.0)), str(dir_dop / f"{stem}_DoP.png"))
            _imwrite(colorjetmap(aop / np.pi), str(dir_aop / f"{stem}_AoP.png"))

            processed_count += 1
            print(f"[Done] image={image_name}, elapsed_sec={time.perf_counter() - t_img:.3f}")

        except KeyboardInterrupt:
            raise
        except Exception as e:
            print(f"[ERROR] image={image_name}, reason={str(e)}")
            traceback.print_exc()
            continue

    elapsed = time.perf_counter() - t0
    print(
        f"summary: total={len(indices)}, processed={processed_count}, skipped={skipped_count}, elapsed_sec={elapsed:.3f}"
    )


if __name__ == "__main__":
    main()
