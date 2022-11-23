function [imgLoG]=filterLoG23(img, sigma, ws)
%Obj: Laplacian of Gaussian filter for 2D or 3D image using convolution
%Input:
%   img: the image
%   sigma: sigma of Gaussian, for 3D [sigma_xy, sigma_z]
%   ws: window size, for 3D [ws_xy, ws_z]

nd=ndims(img);
if ~isa(img,'double');
    img=double(img);
end
if nd==2
    if nargin<3
        ws=ceil(4*sigma);
    end
    x=-ws:ws;
    g = exp(-x.^2/(2*sigma^2));
    % convolutions
    imgXT = padarrayXT(img, [ws ws], 'symmetric');
    fg = conv2(g', g, imgXT, 'valid');  
    % Laplacian of Gaussian
    gx2 = g.*x.^2;
    imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
    imgLoG = imgLoG / (2*pi*sigma^2);
elseif nd==3
    if numel(sigma)==1
        sigma=[sigma sigma];
    end
    if nargin<3
        ws=ceil(2*sigma);
    end
    if numel(ws)==1
        ws=[ws ws];
    end
    % right-hand side of symmetric kernels
    gx = exp(-(0:ws(1)).^2/(2*sigma(1)^2));
    gz = exp(-(0:ws(2)).^2/(2*sigma(2)^2));
    fg = conv3fast(img, gx, gx, gz);
    % Laplacian of Gaussian-filtered input
    gx2 = (0:ws(1)).^2 .*gx;
    gz2 = (0:ws(2)).^2 .*gz;
    fgx2 = conv3fast(img, gx2, gx, gz);
    fgy2 = conv3fast(img, gx, gx2, gz);
    fgz2 = conv3fast(img, gx, gx, gz2);
    imgLoG = (2/sigma(1)^2+1/sigma(2)^2)*fg - ((fgx2+fgy2)/sigma(1)^4 + fgz2/sigma(2)^4);
else
    disp('filterLoG23 only works for 2D and 3D images');
    imgLoG=[];
end
    
end