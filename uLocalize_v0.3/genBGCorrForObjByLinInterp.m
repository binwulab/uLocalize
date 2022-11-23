function [imgCorr,bw2,pixelIdxList2,bndBox,bgCrop,bgPlane_mean,bg_sd]= ...
    genBGCorrForObjByLinInterp(img,bw,pixelIdxList,thickness,bg_ext, bgCorrMethod)
%Objective: generate background correction image around an object
%through linear interpolation of the image around the target point, good
%for both 2D and 3D images
%Input:
%   img: the raw image
%   bw: mask for objects: 1 for all other objects except the current one.
%       the 1 pixel is ignored during bg calculation
%   pixelIdxList: linear indices for the object
%   thickness: the thickness of around the region to be cut, used for bg correction
%   bg_ext: the number of pixel extended from the boundary of the object
%Output: 
%   imgCorr: the background corrected image
%   pixelIdxList2: the linear index of the current obj in the cropped image
%   bndBox: the bounding box for the cropped image [ymin, xmin; ymax, xmax]
%   bgCrop: the linear interpolated plane image around the obj
%   bgPlane_mean: the mean value of interpolated bg plane
%   bw2: the mask file in the cropped image, containing both the current and other objects
%History   
%   B.W., Jan, 2021

numdim=ndims(img);
if nargin<6
    bgCorrMethod = 'plane';
end
[bgData, bndBox, pixelIdxList2, bw2]=getBGDataAroundObj(img,bw,pixelIdxList,thickness,bg_ext);
[par,bg_sd]=fitBGImgWrapper(bgData, bgCorrMethod); %linear fit to the bg image, return para and the std of the difference between the fit and the raw bg
ymin=bndBox(1,1);
ymax=bndBox(2,1);
xmin=bndBox(1,2);
xmax=bndBox(2,2);
if numdim==2
    [X,Y]=meshgrid(xmin:xmax, ymin:ymax);   %gen the coordinate of the cropped image
    bgCrop=genBGImgWrapper(par, Y, X);
    %bgCrop=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3);    %bg correction in the imaging area through interpolation
    imgCorr=double(img(ymin:ymax,xmin:xmax))-bgCrop;    %background corrected image
else
    zmin=bndBox(1,3);
    zmax=bndBox(2,3);
    [X,Y,Z]=meshgrid(xmin:xmax, ymin:ymax,zmin:zmax); %gen the coordinate of the cropped image
    %bgCrop=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3)*(Z-0.5)+par(4); %bg correction in the imaging area through interpolation
    bgCrop=genBGImgWrapper(par, Y, X, Z);    
    imgCorr=double(img(ymin:ymax,xmin:xmax,zmin:zmax))-bgCrop; %background corrected image
end
bgPlane_mean=mean(bgCrop(:));    %calculated mean bg for the interpolated plane
end