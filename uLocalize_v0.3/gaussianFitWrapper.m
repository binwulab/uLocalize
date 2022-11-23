function [pt, status]=gaussianFitWrapper(Im, pt0, fitType, bw)
%The Gaussian Fit Wrapper program for cropped image
%   Imax*exp(-((x-xc)^2+(y-yc)^2)/(2*sigma^2))+bg
%Input
%       Im: the cropped image
%       pt0: the starting point
%       fitType: 'standard gaussian'
%                'integrated gaussian'
%                'integrated gaussian no bg';
%       bw: the mask of other obj: >=1 if the pixel belongs to another object,it will
%           be excluded. The current obj, however, should be set 0. 

%Output
%       pt: [yc, xc, zc, Int, sxy, sz, bg] for 3D
%           [yc, xc, Int, sxy, bg] for 2D
%History
%   BW, Feb 2021
numdim=ndims(Im);
if nargin<3 
    fitType='integrated gaussian';
end
if numdim==3
    [pt, status]=gaussianFit3D(Im, pt0, fitType, bw);
else
    [pt, status]=gaussianFit2D(Im, pt0, fitType, bw);
end
end


