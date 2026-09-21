function [S0,DoLP,AoLP] = cal_Stokes(im)
imOut0 = im(:,:,1:3);
imOut45 = im(:,:,4:6);
imOut90 = im(:,:,7:9);
imOut135 = im(:,:,10:12);

S0 = (imOut0+imOut45+imOut90+imOut135)*0.5;
S1 = imOut0-imOut90;
S2 = imOut45-imOut135;
DoLP = sqrt(S1.^2+S2.^2)./(S0+eps);
AoLP = 0.5*(atan2(S2, S1));
% angles = 0.5*(atan2(S2, S1));
% AoLP = rad2deg(angles);
% AoLP(AoLP < 0) = AoLP(AoLP < 0) + 360;
end