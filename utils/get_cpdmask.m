function mask = get_cpdmask(h,w,pattern)
pr = find((pattern == 'r') + (pattern == 'R'));
num(pr) = 1;
pg = find((pattern == 'g') + (pattern == 'G'));
num(pg) = 2;
pb = find((pattern == 'b') + (pattern == 'B'));
num(pb) = 3;

mask_rgb = zeros(h/2,w/2,3);
mask_rgb(1:2:h/2, 1:2:w/2, num(1)) = 1;
mask_rgb(1:2:h/2, 2:2:w/2, num(2)) = 1;
mask_rgb(2:2:h/2, 1:2:w/2, num(3)) = 1;
mask_rgb(2:2:h/2, 2:2:w/2, num(4)) = 1;

% polarization pattern : [90,45;135,0]
mask = zeros(h,w,12);
for i = 1:3
    mask(2:2:h,2:2:w,i)   = mask_rgb(:,:,i);
    mask(1:2:h,2:2:w,i+3) = mask_rgb(:,:,i);
    mask(1:2:h,1:2:w,i+6) = mask_rgb(:,:,i);
    mask(2:2:h,1:2:w,i+9) = mask_rgb(:,:,i);
end