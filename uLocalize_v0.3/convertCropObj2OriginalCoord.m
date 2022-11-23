function pixelIdxList=convertCropObj2OriginalCoord(pixelIdxList2, bndBox, imSize)
%given an object pixelIdxList in the cropped image given by bndBox,
%calcualte the pixelIdxList in the original image
%Input: 
%   pixelIdxList2: the object in the cropped image
%   bndBox: the bound box of the cropped image
%   imSize: the size of the original image
%Output
%   pixelIdxList: the object in the original image coordinate
%History
%   BW, Jan 2021
%   There is no error checking, assuming the bndBox is within the imSize
numdim=numel(imSize);
ymin=bndBox(1,1);
ymax=bndBox(2,1);
xmin=bndBox(1,2);
xmax=bndBox(2,2);
ny=ymax-ymin+1;
nx=xmax-xmin+1;
if numdim==3
    zmin=bndBox(1,3);
    zmax=bndBox(2,3);
    nz=zmax-zmin+1;
    [yy,xx,zz]=ind2sub([ny,nx,nz],pixelIdxList2);
    yy=yy+ymin-1;
    xx=xx+xmin-1;
    zz=zz+zmin-1;
    pixelIdxList=(zz-1)*imSize(1)*imSize(2)+(xx-1)*imSize(1)+yy;
else
    [yy,xx]=ind2sub([ny,nx],pixelIdxList2);
    yy=yy+ymin-1;
    xx=xx+xmin-1;
    pixelIdxList=(xx-1)*imSize(1)+yy;
end

end