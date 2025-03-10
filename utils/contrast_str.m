% This code is another version of 'graystretch' from "Sparse representation
% -based demosaicing method for microgrid polarimeter imagery"

function Y = contrast_str(X,a)
Y = X;
Y(X>a) = a;
Y(X<0) = 0;
Y = Y./a;
end

% function img=graystretch(img)
% sigma_pd=std(img(:));
% mean_pd=mean(img(:));
% img(img<mean_pd-3*sigma_pd)=mean_pd-3*sigma_pd;
% img(img>mean_pd+3*sigma_pd)=mean_pd+3*sigma_pd;
% maxvalue=max(img(:));
% minvalue=min(img(:));
% img=(img-minvalue)*255/(maxvalue-minvalue);
% end