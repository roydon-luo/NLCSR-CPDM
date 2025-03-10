function Bay_demosaic = bayer_bilinear(mosaic)

conv_kernel_g = [1, 2, 1; 2, 4, 2; 1, 2, 1]/8;
conv_kernel_r = [1, 0, 1; 0, 4, 0; 1, 0, 1]/4;          
           
R_1 = conv2(mosaic(:,:,1), conv_kernel_r, "same");
R = conv2(R_1, conv_kernel_g, "same");
G = conv2(mosaic(:,:,2), conv_kernel_g, "same");
B_1 = conv2(mosaic(:,:,3), conv_kernel_r, "same");
B = conv2(B_1, conv_kernel_g, "same");

Bay_demosaic = cat(3,R,G,B);





