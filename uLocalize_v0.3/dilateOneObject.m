function [pixelIdxList2,pixelIdxDiff]=dilateOneObject(pixelIdxList,imSize,numDilation)
%Objective: dilate objects defined in the connected component structure,
%the pixel value will record the number of objects containing that pixel
%Input
%   pixelIdxList: the pixelIdxList of one object
%   imSize: the image size
%   numDilation: the number of times to gradually increase the size of object
%Output
%   pixelIdxList2: the updated cc structure, each object now contains dilated pixels, in the original coordinate
%   pixelIdxDiff: the difference of pixelIdx of the object before and after dilation in the original coordinate
%History
%   BW: Mar 2021

%% calculate the structure elements for dilation, to ensure adding least number of pixels to the object for each dilation
%making the dist variable persistent. So calculated only once or when
%necessary if dimension of the problem changed
persistent dist numDilation_persist strel;
calc_strel = false;
numdim=numel(imSize);
if isempty(dist) || (numdim==2 && numel(dist)~=50) || (numdim==3 && numel(dist)~=144) %during initiation or num dimension changed
    x=0:9;
    if numdim==2
        [x,y]=meshgrid(x,x);
        dist=sqrt(x.^2+y.^2);    
        clear('x','y');
    else
        [x,y,z]=meshgrid(x,x,x);
        dist=sqrt(x.^2+y.^2+z.^2);
        clear('x','y','z');
    end
    dist=unique(sort(dist,'ascend'));
    dist=dist(2:end);
    calc_strel = true;
end
if numDilation > numel(dist)
    disp('dilateObjects: you have tried to dilate to much, dear');
    return;
end
radius=dist(numDilation);
radiusInt=floor(radius);

%%save some time not to calculate this strel every time. 
if isempty(numDilation_persist)
    numDilation_persist=numDilation;
    calc_strel=true;
end
if numDilation_persist ~= numDilation
    calc_strel = true;
    numDilation_persist = numDilation;
end

if calc_strel
    % calculate the structure element for dilation
    if numdim==2
        [xx,yy]=meshgrid(1:2*radiusInt+1, 1:2*radiusInt+1);
        strel=sqrt((xx-radiusInt-1).^2+(yy-radiusInt-1).^2)<radius;
        clear('xx','yy');
    else
        [xx,yy,zz]=meshgrid(1:2*radiusInt+1, 1:2*radiusInt+1, 1:2*radiusInt+1);
        strel=sqrt((xx-radiusInt-1).^2+(yy-radiusInt-1).^2+(zz-radiusInt-1).^2)<=radius;
        clear('xx','yy','zz');
    end
end
%% dilate each object, with possible objects overlapping in pixel, but each object is still separate

%bw2=bw;     %This is too inefficient to do object by object, needs to improve later
radiusInt=ceil(radius);
    %generate cropped binary image just containing the current object
    %Note, I have not given the bw in the getObjBndBox, so bwTmp will be
    %binary image containing only the current object
    [bndBox,~,bwTmp]=getObjBndBox(pixelIdxList, imSize, radiusInt);
    bwTmp2=imdilate(bwTmp,strel);    %dilate the cropped binary image
    pixelIdxList2=find(bwTmp2);
    pixelIdxDiff=find(bwTmp2 & ~bwTmp); %find the expanded pixel
    pixelIdxList2=convertCropObj2OriginalCoord(pixelIdxList2,bndBox,imSize);
    pixelIdxDiff=convertCropObj2OriginalCoord(pixelIdxDiff,bndBox,imSize);
    %Note: each pixel value will be a count of how many object containing it
    %bw2(pixelIdxDiff)=bw2(pixelIdxDiff)+1;    
    %Do we need to update the Object characteristic here since the obj is dilated? 
% if cc.nObjects(1)~=0
%     areaDiff=mean(areaDiff(1:cc.nObjects(1)));  %average the singlet
% else
%     areaDiff=mean(areaDiff);
% end
end
