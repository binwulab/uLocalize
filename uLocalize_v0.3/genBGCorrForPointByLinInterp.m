function [imgCorr, fgImg2, pt2, bndBox, bgCrop, bgPlane_mean, bg_sd, pixelIdxList2, bg_ext]= ...
    genBGCorrForPointByLinInterp(img, fgImg, pt, cutwidth, thickness, ROISizeOption,bgCorrMethod, pixelIdxList)
%Objective: generate background correction image around a point
%through linear interpolation of the image around the target point, good
%for both 2D and 3D images. The code is used for both LM and CC mode
%Input:
%   img: the raw image
%   fgImg: foreground image, count the number of obj a pixel belongs to
%   pt: the center of the current object: [Y,X,Z]
%   cutwidth: the size of the image: [cutwidth_xy, cutwidth_z]
%           image returned: [2*cutwidth+1, 2*cutwidth+1]
%   thickness: the thickness of around the region to be cut, used for bg correction
%   Optional: 
%   ROISizeOption: 'small': only the center pixels around the object, default
%       'large': include the surrounding pixels,
%   bgCorrMethod: background correction method, 'plane', 'quad', 'median'
%   pixelIdxList: the linear pixel Index of the currentd object in the
%       original image coordinates, only used for the CC mode
%Output: 
%   imgCorr: the background corrected image
%   fgImg2: the cropped foregroung image, may contain both the current and other objects
%   pt2: the point transformed to the croppedcoordinate
%   bndBox: the bounding box for the cropped image [ymin, xmin; ymax, xmax]
%   bgCrop: the linear interpolated plane image around the obj
%   pixelIdxList2: the linear index of the current obj in the cropped image
%   bgPlane_mean: the mean value of interpolated bg plane
%   bg_sd: the std of difference between bg image and the plane fit 
%   bg_ext: the distance between the boundary of the object and the ROI, used for background correction
%History   
%   B.W., Jan, 2021
%   B.W., Aug, 2021
    %Change the code for converting pixelIdx to the cropped coord, the bw
    %   matrix will not change, therefore should be more memory efficient
    %Unify this code for both LM and CC modes. For LM mode, note that the
    %   fgImg2 returned containing the number of object


if nargin<5
    disp('Error, genBGCorrForPointByLinInterp: expect at least 5 arguments');
    return;
end

if nargin<6
    ROISizeOption = 'small';
end

if nargin<7
    bgCorrMethod = 'plane';
end

sz=size(img);
numdim=ndims(img);
data=getBGDataAroundPoint(img,fgImg,pt,cutwidth,thickness);
[par, bg_sd]=fitBGImgWrapper(data,bgCorrMethod);
bndBox=getPointBndBox(pt, cutwidth, thickness, sz, ROISizeOption); %generate the bounding box around the curr obj for actual fitting
ymin=bndBox(1,1);
ymax=bndBox(2,1);
xmin=bndBox(1,2);
xmax=bndBox(2,2);
yc2=ceil(pt(1))-ymin+0.5;
xc2=ceil(pt(2))-xmin+0.5;
if numdim == 2
    pt2=[yc2,xc2];
    [X,Y]=meshgrid(xmin:xmax, ymin:ymax);
    %bgCrop=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3);    %bg correction in the imaging area through interpolation
    bgCrop=genBGImgWrapper(par, Y, X);
    imgCorr=double(img(ymin:ymax,xmin:xmax))-bgCrop;    %background corrected image
    %The fg counts around current pt is reduced by 1 throughout the cropped
    %box. Note fgImg is uint8. so at most 1 will be reduced to 0, 0 still will be 0.
    %fgImg2=fgImg(ymin:ymax,xmin:xmax)-1;    
    fgImg2=fgImg(ymin:ymax,xmin:xmax);    % the cropped foreground image
else
    zmin=bndBox(1,3);
    zmax=bndBox(2,3);
    zc2=ceil(pt(3))-zmin+0.5;
    pt2=[yc2,xc2,zc2];
    [X,Y,Z]=meshgrid(xmin:xmax, ymin:ymax, zmin:zmax);
    %bgCrop=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3)*(Z-0.5)+par(4);
    bgCrop=genBGImgWrapper(par, Y, X, Z);    
    imgCorr=double(img(ymin:ymax,xmin:xmax,zmin:zmax))-bgCrop; %background corrected image
    %fgImg2=fgImg(ymin:ymax,xmin:xmax,zmin:zmax)-1;
    fgImg2=fgImg(ymin:ymax,xmin:xmax,zmin:zmax);    % the cropped foreground image
end
bgPlane_mean=mean(bgCrop(:));    %calculated mean bg for the interpolated plane

% Now convert pixelInxList from the original image coord to the cropped image
% Note I did not change the fgImg at all, which is important to maintain
% low memory consumption
if nargout>7 
    if ~isempty(pixelIdxList)  %if the object pixelIdxList is given
        if numdim==2
            [yy,xx]=ind2sub(sz, pixelIdxList);  %convert linear pixel index to coordinate
            yy=yy-ymin+1;   %convert to the cropped coordinate, a list of y coordinates for the pixels of the object
            xx=xx-xmin+1;
            ny=ymax-ymin+1;
            nx=xmax-xmin+1;
            valid=(yy>=1 & yy<=ny & xx>=1 & xx<=nx);   %Important: keep only the points in the bounding box
            yy=yy(valid);
            xx=xx(valid);
            pixelIdxList2=(xx-1)*ny+yy;
            %bg_ext is the distance between the object and the boundary of the ROI
            bg_ext=mean([min(yy),min(xx),ny-max(yy),nx-max(xx)]);
        else
            [yy,xx,zz]=ind2sub(sz, pixelIdxList);   %convert linear pixel index to coordinate
            yy=yy-ymin+1;   %convert to the cropped coordinate, a list of y coordinates for the pixels of the object
            xx=xx-xmin+1;
            zz=zz-zmin+1;
            ny=ymax-ymin+1;
            nx=xmax-xmin+1;
            nz=zmax-zmin+1;
            valid=(yy>=1 & yy<=ny & xx>=1 & xx<=nx & zz>=1 & zz<=nz);   %Important: keep only the points in the bounding box
            yy=yy(valid);
            xx=xx(valid);
            zz=zz(valid);
            pixelIdxList2=(zz-1)*ny*nx+(xx-1)*ny+yy;
            %bg_ext is the distance between the object and the boundary of the ROI
            bg_ext=zeros(1,2);
            bg_ext(1)=mean([min(yy),min(xx),ny-max(yy),nx-max(xx)]);
            bg_ext(2)=mean([min(zz), nz-max(zz)]);
        end
    else
        pixelIdxList2=[];
        bg_ext=[];
    end
end