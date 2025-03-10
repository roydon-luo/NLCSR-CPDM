# Learning a Non-Locally Regularized Covolutional Sparse Representation for Joint Chromatic and Polarimetric Demosaicking
(2025-3-10) Our codes are released!
## 1.Testing
For demosaicking of real images captured by DoFP camera, please put mosaic data in:
```
'data/real/'
```
and run:
```bash
demo_real.m
```
For demosaicking of ground truth captured by DoTP camera, please put full-resolution data (0°,45°,90°,135°) in:
```
'data/synthetic/your_selected_image/'
```
and run:
```bash
demo_syn.m
```
## 2.Parameters
The setting of parameters affects the final outcome of demosaicking, and you can adjust them according to your own needs. For instance, when the noise level of DoLP is high, the parameter c, which influences the non-local self-similarity constraint, should be correspondingly increased. Conversely, it should be decreased to avoid mistakenly removing details as artifacts.
##
Welcome to adjust the structure and parameters for improvement of NLCSR-CPDM's performance!
## 3.Citation
If you find our work useful and use it for your research, please cite our paper:
```
@articles{luo2024learning,
  title={Learning a Non-Locally Regularized Covolutional Sparse Representation for Joint Chromatic and Polarimetric Demosaicking},
  author={Luo, Yidong and Zhang, Junchao and Shao, Jianbo and Tian, Jiandong and Ma, Jiayi},
  journal={IEEE Transactions on Image Processing},
  volume={33},
  pages={5029--5044},
  year={2024}
}
```

