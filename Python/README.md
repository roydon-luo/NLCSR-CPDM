# Python implementation

This directory is a standalone Python implementation of NLCSR-CPDM.  It supports CPU execution by default and optional CUDA acceleration through CuPy.  It does **not** require MATLAB.

## Requirements

- Python 3.9 or later
- NumPy, SciPy, and ImageIO
- Optional: CuPy matching your CUDA version for GPU acceleration

Install the required CPU packages from this directory:

```bash
python -m pip install -r requirements.txt
```

For CUDA, install the CuPy wheel appropriate for your local CUDA toolkit; see the [CuPy installation guide](https://docs.cupy.dev/en/stable/install.html).  Omitting `--use-cuda` always uses the CPU implementation.

## Quick start: real DoFP demosaicking

Run from this `Python` directory.  The included sample images are in `data/real/`.

```bash
python demo_real_save.py --image-root data/real --save-dir real_results
```

To use CUDA on GPU 0:

```bash
python demo_real_save.py --image-root data/real --save-dir real_results --use-cuda --gpu-id 0
```

The program recursively accepts PNG input under `--image-root` and saves reconstructed polarization channels, Stokes-derived images, and MATLAB-compatible `.mat` output under `--save-dir`.

## Grayscale polarization denoising

The denoising example accepts common image formats and recursively scans its input directory.  A small example is included in `data/Noise/Owl/`.

```bash
python demo_real_Denoise.py --image-root data/Noise/Owl --save-dir denoise_results
```

Use `--use-cuda --gpu-id 0` to enable CUDA.

## Multi-GPU processing

Both demos have multi-GPU launchers.  They split an input tree across the specified GPUs and merge the outputs when complete.

```bash
python run_demo_real_save_multigpu.py --gpus 0,1 --image-root data/real --save-dir real_results_multi --use-cuda
python run_demo_real_Denoise_multigpu.py --gpus 0,1 --image-root data/Noise/Owl --save-dir denoise_results_multi --use-cuda
```

## Input layout

- `data/real/`: DoFP mosaic images for `demo_real_save.py`.
- `data/Noise/`: grayscale polarization images for `demo_real_Denoise.py`.
- `data/synthetic/`: four-angle reference sets retained as example data.

Generated output directories, Python bytecode, and CuPy caches are intentionally ignored by Git.  Algorithm parameters are defined in `para_set.py`; adjust them only after establishing a baseline run.

## Citation

If you use the Python implementation or the related dataset in your research, please cite:

```bibtex
@inproceedings{song2026pol,
  title={Pol-CACTI: A System and Dataset for High-Speed Polarized Video Compressive Imaging},
  author={Song, Yunfeng and Luo, Yidong and Wang, Ping and Jiang, Xingjian and Yuan, Xin},
  booktitle={European Conference on Computer Vision},
  pages={521--539},
  year={2026},
  organization={Springer}
}
```
