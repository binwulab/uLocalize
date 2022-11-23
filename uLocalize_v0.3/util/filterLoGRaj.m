function imgLoG=filterLoGRaj(img,sigma_xy,ws)
%Obj: The Laplacian of Gaussian filter from Raj
%Input: 
%   img: the raw image
%   sigma_xy, the sigma_xy for the LoG filter
%   ws: window size, the size of the filter

if nargin<3
    ws=max(ceil(4*sigma_xy+1),5);   %at least 5 pixels filter size
    if mod(ws,2)==0     %make it an odd number
        ws=ws+1;
    end
end
nd=ndims(img);
op=fspecial('log',ws,sigma_xy);
op=op-sum(op(:))/numel(op);
if nd==3
    op=1/3*cat(3,op,op,op); %actually just concatenate 3 2D LoG filter in 3D
end
imgLoG=-imfilter(double(img),op,'symmetric');
end