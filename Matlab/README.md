# MATLAB implementation

This is the original MATLAB implementation of NLCSR-CPDM.  Run MATLAB with the current folder set to this `Matlab` directory; the demos add the required `utils` directories automatically.

## Requirements

- MATLAB with Image Processing Toolbox (`imread`, `imshow`, and `im2double` are used)

## Real DoFP demosaicking

Place DoFP mosaic PNG files in `data/real/` and run:

```matlab
demo_real
```

The included `album.png` and `carlight.png` files can be used directly.  Results are written to `real_results/`.

## Synthetic DoTP evaluation

Each reference scene must have its full-resolution polarization images at:

```text
data/synthetic/<dataset>/<scene>/0.png
data/synthetic/<dataset>/<scene>/45.png
data/synthetic/<dataset>/<scene>/90.png
data/synthetic/<dataset>/<scene>/135.png
```

The included `Monno`, `Qiu`, and `Wen` datasets follow this layout.  Select the dataset in `demo_syn.m` (`dataset_name`) and run:

```matlab
demo_syn
```

The synthetic demo uses a cropped region for quick testing.  Remove or change `S_ori = S_ori(301:516,401:616,:);` in `demo_syn.m` to process another region or the full image.  Results are written to `synthetic_results/`.

## Parameters

Algorithm parameters are defined in `para_set.m`.  In particular, `C` controls the non-local self-similarity constraint: increasing it can suppress high DoLP noise, while decreasing it can help preserve fine details.  Establish a baseline result before changing parameters.

## Citation

If you find this work useful, please cite:
```
@articles{luo2024learning,
  title={Learning a Non-Locally Regularized Covolutional Sparse Representation for Joint Chromatic and Polarimetric Demosaicking},
  author={Luo, Yidong and Zhang, Junchao and Shao, Jianbo and Tian, Jiandong and Ma, Jiayi},
  journal={IEEE Transactions on Image Processing},
  volume={33},
  pages={5029--5044},
  year={2024},
  publisher={IEEE}
}
```
