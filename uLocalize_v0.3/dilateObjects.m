function [fgImg,cc2,areaDiff]=dilateObjects(bw,cc,numDilation)
%Objective: dilate objects defined in the connected component structure,
%the pixel value will record the number of objects containing that pixel
%Input
%   bw: the beginning mask file
%   cc: the connected component strucutre
%   numDilation: the number of times to gradually increase the size of object
%Output
%   fgImg: the updated mask file, with each pixel recording the number of obj containing that pixel
%   cc2: the updated cc structure, each object now contains dilated pixels
%   areaDiff: [1,NumObjects], the difference of pixel number in each object before and after dilation
%History
%   BW: Jan 2021
%   BW: Aug 2021, 
        %change BW2 from double to uint8, which increases speed and reduce memory significantly
        %try to make strel as persistent variable, probably not a good idea for multi-thread execution and not saving much time

%% calculate the structure elements for dilation, to ensure adding least number of pixels to the object for each dilation
%making the dist variable persistent. So calculated only once or when
%necessary if dimension of the problem changed
persistent dist numDilation_persist strel;
calc_strel = false;
numdim=ndims(bw);
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
cc2=cc;
fgImg=zeros(size(bw),'uint8');
radiusInt=ceil(radius);
areaDiff=zeros(1,cc.NumObjects);
for i=1:cc.NumObjects   %Dilate all objects
    %generate cropped binary image just containing the current object
    %Note, I have not given the bw in the getObjBndBox, so bwTmp will be
    %binary image containing only the current object
    [bndBox,~,bwTmp]=getObjBndBox(cc.PixelIdxList{i}, size(bw), radiusInt);
    bwTmp=imdilate(bwTmp,strel);    %dilate the cropped binary image
    pixelIdxList2=find(bwTmp);
    cc2.PixelIdxList{i}=convertCropObj2OriginalCoord(pixelIdxList2,bndBox,size(bw));
    %Note: each pixel value will be a count of how many object containing it
    fgImg(cc2.PixelIdxList{i})=fgImg(cc2.PixelIdxList{i})+1;    
    areaDiff(i)=numel(cc2.PixelIdxList{i})-numel(cc.PixelIdxList{i});
    %Do we need to update the Object characteristic here since the obj is dilated? 
end
% if cc.nObjects(1)~=0
%     areaDiff=mean(areaDiff(1:cc.nObjects(1)));  %average the singlet
% else
%     areaDiff=mean(areaDiff);
% end
end
