function [ISum, dist2ObjCtr, dist2Obj]=calcIntSumGranule(img, p, cc, ts)
%Objective: given a granule structure cc, calculate the intensity of a
%given image for these granules
%Input:
%   img: the raw image, typically a different channel from the detected cc
%   p: the detection parameter structure
%   cc: the cc structure for the granules
%   ts: the matrix containing the fitted ts: [y,x,z,int,sigma_xy,sigma_z],
%       the number of rows should equal the number of obj in cc
%Output
%   ISum: the ISum of the granules in the current image
%   dist2ObjCtr: distance of the ts to the object center
%   dist2Obj, distance of the ts to the object boundary
%History:
% B.W. Mar, 2021
%      Apr, 2021: Add input para ts and output paras dist

if ~isequal(size(img), cc.ImageSize) || cc.NumObjects ~= size(ts,1)
    disp('calcIntSumGranule: image size does not match with cc');
    ISum=[];
    return
end
if nargin>3     %calculate the distance only when the ts is given
    if nargout>1    
        dist2ObjCtr=zeros(cc.NumObjects,1);
    end
    if nargout>2 
        dist2Obj=zeros(cc.NumObjects,1);
    end
end
numdim=ndims(img);
bw=zeros(size(img));
%set all thresholded pixels to 1. At this moment, the object should not
%overlap. After background extension, it may. We haven't taken that into
%account. 
PixIdxList=cell2mat(cc.PixelIdxList');
bw(PixIdxList)=true;    

%Now taken into account that the objects might overlap
% for i=1:cc.NumObjects
%     bw(cc.PixelIdxList{i})=bw(cc.PixelIdxList{i})+1;
% end

ISum=zeros(cc.NumObjects,1);


thickness=p.thickness;
bg_ext=p.bg_extension;

%%calculate ISum
for i=1:cc.NumObjects
    pixelIdxList=cc.PixelIdxList{i};
    [imgCorr, bw2,pixelIdxList2, bndBox, bgCrop, bgPlane_mean, bgCorr_sd]= ...
        genBGCorrForObjByLinInterp(img,bw,pixelIdxList,thickness,bg_ext, 'quad');
    ISum(i)=sum(imgCorr(pixelIdxList2)./double(bw2(pixelIdxList2))); %if one pixel belongs 2 objects, divide the intensity by 2
%    ISum(i)=sum(imgCorr(pixelIdxList2)); %This ignores the object overlap
end

%%Calculate the distance to the ObjCtr
if nargout>1
    for i=1:cc.NumObjects
        dist2ObjCtr(i)=sqrt(sum((cc.ObjCtr(i,1:numdim)-ts(i,1:numdim)).^2));
    end
end

%%Calculate the distance to the Obj
if nargout>2
    imSize=cc.ImageSize;
    for i=1:cc.NumObjects
        pixelIdxList=cc.PixelIdxList{i};
        tsPos=ts(i,1:numdim);
        if numdim == 3
            [yobj,xobj,zobj] = ind2sub(imSize,pixelIdxList);
            objPos=[yobj,xobj,zobj];
        else
            [yobj,xobj] = ind2sub(imSize,pixelIdxList);
            objPos=[yobj,xobj];
        end
        tsPos=repmat(tsPos, numel(pixelIdxList),1);
        dist2Obj(i)=sqrt(min(sum((objPos-tsPos).^2,2)));
    end
end
end