function pixelIdxList2=convertCoordOriginal2Crop(pixelIdxList, bndBox, imSize)
%given an object in an image, calculate the pixelIdxlist in the cropped
%image given by the bndBox
%Input: 
%   pixelIdxList: the object in the original image
%   bndBox: the bnd box of the cropped image 
%   imSize: the size fo the original image
%Output:
%   pixelIdxList2: the object in the cropped image coordinate
%History
%   BW, Aug, 2021
numdim=numel(imSize);
ymin=bndBox(1,1);
ymax=bndBox(2,1);
xmin=bndBox(1,2);
xmax=bndBox(2,2);
ny=ymax-ymin+1;
nx=xmax-xmin+1;
if numdim ==3 
    zmin=bndBox(1,3);
    zmax=bndBox(2,3);
    nz=zmax-zmin+1;
    [YY,XX,ZZ]=ind2sub(imSize, pixelIdxList); %coord in the original image
    yy=YY-ymin+1;
    xx=XX-xmin+1;
    zz=ZZ-zmin+1;
    pixelIdxList2=(zz-1)*ny*nz+(xx-1)*ny+yy;
else
    [YY,XX]=ind2sub(imSize, pixelIdxList); %coord in the original image
    yy=YY-ymin+1;
    xx=XX-xmin+1;
    pixelIdxList2=(xx-1)*ny+yy;
end
end