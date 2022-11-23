function [final_pts, p, detection_res,dummy]=uLocalizeFitGranule(img, p, cc, verbose)
%Obj: Fit the objects detected by connected component of thresholded image
%Input:
%   img: the smoothed image
%   p: the detection parameter structure
%   cc: the connected component using thresholded predetection
%   verbose: output the intermediate information
%Output
%   final_pts: [y,x,c,Int] of the granule, the center and intensity
%   p: the detection parameter structure
%   detection_res: the cc structure after dilating the image
%History
%   BW: Mar 2021

if verbose
    disp('uLocalizeFitGranule: beginning ...');
    tic;
end
numdim=ndims(img);
bw=false(size(img));
%set all thresholded pixels to 1. At this moment, the object should not
%overlap. After background extension, it may. 
PixIdxList=cell2mat(cc.PixelIdxList');
bw(PixIdxList)=true;    
% bwCount=bw;
cc.ISum=zeros(cc.NumObjects,1);
thickness=p.thickness;
bg_ext=p.bg_extension;
for i=1:cc.NumObjects
    pixelIdxList0=cc.PixelIdxList{i};
    %generate a cropped, local background corrected image
    [imgCorr, ~,pixelIdxList2, bndBox, bgCrop, bgPlane_mean, bgCorr_sd]= ...
        genBGCorrForObjByLinInterp(img,bw,pixelIdxList0,thickness,bg_ext, p.bgCorrMethod);
%     INormByObjNum=imgCorr(pixelIdxList2)./ double(bw2(pixelIdxList2));
%     ISum=sum(INormByObjNum(:));
    ISum=sum(imgCorr(pixelIdxList2));
    numPix=numel(pixelIdxList2);
    growGranule=true;
    numDilation=1;
    while growGranule && ((numdim==3 && numDilation<144) || (numdim==2 && numDilation<50))
        [pixelIdxList1, pixelIdxDiff1]=dilateOneObject(pixelIdxList0,size(bw), numDilation);
        numPix1=numel(pixelIdxList1);
        [imgCorr, bw2,pixelIdxList2, bndBox, bgCrop, bgPlane_mean, bgCorr_sd]= ...
            genBGCorrForObjByLinInterp(img,bw,pixelIdxList1,thickness,bg_ext, p.bgCorrMethod);
        %bwCount(pixelIdxDiff)=bw(pixelIdxDiff)+1;
%         INormByObjNum=imgCorr(pixelIdxList2)./double(bw2(pixelIdxList2));   %if one pixel belongs 2 objects, divide the intensity by 2
%         ISum2=sum(INormByObjNum(:));          %somehow, this does not converge, so I did not use this
        ISum2=sum(imgCorr(pixelIdxList2)); %This does not consider the overlap between objects
        dI=ISum2-ISum;
        if dI<bgCorr_sd*sqrt(numPix1-numPix)
            growGranule=false;
            break;
        end
        ISum=ISum2;
        numPix=numPix1;
        numDilation=numDilation+1;
    end
    if numDilation>1
        %In this approach, ignore the last round of dilation entirely. 
%         [pixelIdxList1, pixelIdxDiff1]=dilateOneObject(pixelIdxList0,size(bw), numDilation-1);  
%         cc.PixelIdxList{i}=pixelIdxList1;
        %In the 2nd approach, take those with pixel value largr than
        %background+std in the last round of dilation
%         pixelIdxDiff=pixelIdxDiff1(imgCorr(pixelIdxDiff1)>bgCorr_sd);
%         pixelIdxDiff=convertCropObj2OriginalCoord(pixelIdxDiff, bndBox, size(img));
%         cc.PixelIdxList{i}=union(pixelIdxList1,pixelIdxDiff);
        %In the 3rd option, just use numDilation
        ISum=ISum2;
        cc.PixelIdxList{i}=pixelIdxDiff1;
    end
    bw(pixelIdxDiff1)=bw(pixelIdxDiff1)+1;  %This is to handle the case where granule are close to each other
    cc.ISum(i)=ISum;    %Note: the last dilation is ignored
end

%prepare for the output
ObjCtr=regionprops(cc,'Centroid');
ObjCtr = reshape([ObjCtr.Centroid],numdim, cc.NumObjects)';
ObjCtr(:,1:2) = [ObjCtr(:,2),ObjCtr(:,1)];  %change the order of x,y
cc.ObjCtr=ObjCtr;
final_pts=[ObjCtr-0.5, cc.ISum];    %The -0.5 is used to offset the pixel indexing
if verbose
    disp(['uLocalizeFitGranule: Time spent: ', num2str(toc)]);
end
detection_res=cc;
dummy=[];
end

