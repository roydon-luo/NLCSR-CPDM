function [I_jet] = colorjetmap(I)
I_gray = rgb2gray(I);                
I_jet = ind2rgb(gray2ind(I_gray, 256), jet(256));
I_jet = im2uint8(I_jet);
