function [final_pts, p, detection_res, fgImg]=uLocalizeFitCC(img, p, cc, calibrateDilation, verbose)
%Obj: Fit the objects detected by connected component of thresholded image
%Input:
%   img: the raw image
%   p: the detection parameter structure
%   cc: the connected component using thresholded predetection
%   calibration: logical variable for calibration or not
%   verbose: output the intermediate information
%Output
%   final_pts: 
%       for GaussianMask2D: [yc,xc,N0,Chi2]
%       for GaussianFit2D: [yc,xc,N0,sigma_xy]
%       for GaussianMask3D: [yc,xc,zc,N0,err]
%       for GaussianFit3D: [yc,xc,zc,N0,sigma_xy,sigma_z]
%   p: the output detection parameter structure
%   detection_res=struct('single',spotSingle,'multiple',spotMultiple,'cc',cc);
%       single: the single spot detection results
%           fitRes: [y,x,z,Int, err0] for GaussianMask
%           INorm: the normalized intensity using the median of good single spot
%           fitStatus: the fit status parameter for each spot
%               0: N/A
%               1: good fit
%               2: Not converge
%               3: Results falling out of the object
%           ISum: the integratd intensity
%           bg_mean: the mean background for each spot
%           bg_sd: the std of the background after correction with plane image
%           nPix: number of pixel for this object
%           bg_ext: the background extension for each object, i.e. the distance between the object boundary and the cropped image
%       multiple
%           fitRes: [y,x,z,Int]
%           INorm: the normalized intensity using the median of good single spot
%           ISum: the integratd intensity
%           bg_mean: the mean background for each spot
%           bg_sd: the std of the background after correction with plane image
%           nPix: number of pixel for this object
%           nRNA: the number of RNA in each object
%           cc: the cc structure
%   fgImg: the image: each pixel containing the number of object in that pixel
%History
%   BW: Jan 2021
%   BW: Aug 2021: 
        %separate BG spots from multiple, use single fitting to treat them
        %For dilated image, use the type uint8, which is much more efficient
        %Add parfor in uLocalizeFitCCSingle, more efficient for calibration
if nargin<5
    verbose=true;
end
if verbose
    disp('uLocalizeFitCC: beginning ...');
    tic;
end
numdim=ndims(img);
bw=false(size(img));
%set all thresholded pixels to 1. At this moment, the object should not
%overlap. After background extension, it may. 
bw(cell2mat(cc.PixelIdxList'))=true;    

if calibrateDilation    %calibrate the number of dilation needed and bg_extension
    if verbose
        disp('uLocalizeFitCC: calibrating parameters for CC mode');
    end
    p=uLocalizeCalibrateWithSingle(img, p, cc, bw);
end
numDilation=p.numDilation;
bg_extension=p.bg_extension;
ISum2IFitRatio=p.ISum2IFitRatio;
if numDilation >0   %Dilate all objects, update the cc and bw
    [fgImg, cc]=dilateObjects(bw,cc, numDilation); %Note: fgImg is uint8 after calling this line, otherwise, it is logical
    %p.aMax=p.aMax+areaDiff;
else
    fgImg=uint8(bw);
end
clear('bw');

%% fitting single spots
%   210817: discarded spots with fitStatus~=1. So the number of single
%   spots in the spotSingle is different from nSingle r performance for parallelization
nSingle=cc.nObjects(1);
if nSingle ~= 0  %singles
    ObjCtr=cc.ObjCtr(1:nSingle,:);
    pixelIdxList=cc.PixelIdxList(1:nSingle);   %This results in better performance for parallel
    [fitRes, fitStatus, ISum, bg_mean, bg_sd, nPix, bg_ext] = ...
        uLocalizeFitSingleWrapper(img, fgImg, p, ceil(ObjCtr),pixelIdxList);
   
    indValid=find(fitStatus==1);
    ISingle=fitRes(indValid,numdim+1);  %Only use those with good fit to calculate single statistics
    ISingleMed=median(ISingle);   %median of single spot intensity
    %ISingleQuartile=quantile(ISingle,0.25);
    %fitRes(fitStatus>1,numdim+1)=ISum(fitStatus>1);  %if gaussian mask does not converg or out of the object, replace it with integrated intensity
    %INorm=fitRes(:,numdim+1)/ISingleMed;   %Normalize the intensity by the median single intensity
    INorm=ISingle/ISingleMed;
    ns=numel(indValid);
    spotSingle=struct('fitRes',fitRes(indValid,:),'INorm',INorm, 'fitStatus', fitStatus(indValid),...
        'ISum',ISum(indValid),'bg_mean',bg_mean(indValid),'bg_sd',bg_sd(indValid),'nPix',nPix(indValid), 'bg_ext',bg_ext(indValid,:));
else
    spotSingle=struct('fitRes',[],'INorm',[], 'fitStatus', [],...
        'ISum',[],'bg_mean',[],'bg_sd',[],'nPix',[], 'bg_ext',[]);
    %spotSingle=[];
    ns=0;
end
%% deal with the complex spots, only calculate the integrated intensities
nMultiple=cc.nObjects(2);
if nMultiple>0  %210815: change this to only multiples
    ISum=zeros(nMultiple,1);
    bg_mean=zeros(nMultiple,1);
    bg_sd=zeros(nMultiple,1);
    nPix=zeros(nMultiple,1);
    fitRes=zeros(nMultiple,numdim+1);

    for i=1 : nMultiple
        pixelIdxList=cc.PixelIdxList{nSingle+i};
        objCtr=ceil(cc.ObjCtr(nSingle+i,:));
        %generate a cropped, local background corrected image
        [imgCorr, fgImg2,pixelIdxList2, bndBox, bgCrop, bgPlane_mean, bgCorr_sd]= ...
            genBGCorrForObjByLinInterp(img,fgImg,pixelIdxList,p.thickness,bg_extension, p.bgCorrMethod); 
        %calculate the intensity sum, Note in bw2, each pixel value represents how many obj contained that pixel
        INormByObjNum = imgCorr(pixelIdxList2)./double(fgImg2(pixelIdxList2));
        ISum(i)=sum(INormByObjNum(:));
        bg_mean(i)=bgPlane_mean;
        bg_sd(i)=bgCorr_sd;
        nPix(i)=numel(pixelIdxList);
        fitRes(i,1:numdim)=objCtr-0.5;  %Center of object, -0.5 is used to offset the indexing
        fitRes(i,numdim+1)=ISum(i); %The intensity is currently set to the sum of pixel
    end
    intensity=fitRes(:,numdim+1)/ISum2IFitRatio;    %convert integrated intensity to fitted intensity
    INorm=intensity/ISingleMed;
    %nRNA = max([(INorm> 1-ISingleQuartile/ISingleMed & INorm<1.5),round(INorm)],[], 2); %calculate the number of mRNA
    nRNA = max([(INorm> p.ISingleLowBound & INorm<1.5),round(INorm)],[], 2); %calculate the number of mRNA
    indValid=find(nRNA>=1);
    nm=numel(indValid);
    spotMultiple=struct('fitRes',fitRes(indValid,:),'INorm', INorm(indValid), 'ISum', ISum(indValid), ...
        'bg_mean',bg_mean(indValid),'bg_sd',bg_sd(indValid), 'nPix', nPix(indValid),'nRNA',nRNA(indValid));
else
    spotMultiple=struct('fitRes',[],'INorm', [], 'ISum', [], ...
        'bg_mean',[],'bg_sd',[], 'nPix', [],'nRNA',[]);
    nm=0;
end

%% deal with BG spots, which have Area less than aMin but larger than 1. treated it the same way as single
nBG=cc.nObjects(3);
if nBG>0
    ObjCtr=cc.ObjCtr(nSingle+nMultiple+1:end,:);        %This results in bette parallelization
    pixelIdxList=cc.PixelIdxList(nSingle+nMultiple+1:end);   %Note this is a cell array
    [fitRes, fitStatus, ISum, bg_mean, bg_sd, nPix, bg_ext] = ...
        uLocalizeFitSingleWrapper(img, fgImg, p, ceil(ObjCtr),pixelIdxList);
    INorm=fitRes(:,numdim+1)/ISingleMed;
    %nRNA=max([(INorm > 1-ISingleQuartile/ISingleMed & INorm<1.5),round(INorm)],[], 2); %calculate the number of mRNA
    nRNA=max([(INorm > p.ISingleLowBound & INorm<1.5),round(INorm)],[], 2); %calculate the number of mRNA
    indValid=find(nRNA>=1 & fitStatus == 1);
    nb=numel(indValid);
    spotBG=struct('fitRes',fitRes(indValid,:),'INorm',INorm(indValid), 'fitStatus', fitStatus(indValid), ...
        'ISum',ISum(indValid),'bg_mean',bg_mean(indValid),'bg_sd',bg_sd(indValid),'nPix',nPix(indValid), 'bg_ext',bg_ext(indValid,:), 'nRNA', nRNA(indValid));
else
    spotBG=struct('fitRes',[],'INorm',[], 'fitStatus', [], ...
        'ISum',[],'bg_mean',[],'bg_sd',[],'nPix',[], 'bg_ext',[], 'nRNA', []);
    nb=0;
end

%% Make the final_pts
if ns+nm+nb>0
    switch upper(p.method)
        case 'GAUSSIANFIT'  %To be implemented
            final_pts=zeros(ns+nm+nb, 2*numdim); %[yc,xc,N0,sigma_xy]
            if ns>0
                final_pts(1:ns,:)=spotSingle.fitRes;
            end
            if nm>0
                final_pts(ns+1:ns+nm,1:(numdim+1))=spotMultiple.fitRes(1:end,1:(numdim+1));
            end
            if nb>0
                final_pts(ns+nm+1:end,:)=spotBG.fitRes;
            end
        otherwise   %The default is GaussianMask
            final_pts=zeros(ns+nm+nb, numdim+2); %[yc,xc,zc,N0,err]
            if ns>0
                final_pts(1:ns,:)=spotSingle.fitRes;
            end
            if nm>0
                final_pts(ns+1:ns+nm,1:(numdim+1))=spotMultiple.fitRes(1:end,1:(numdim+1));
            end
            if nb>0
                final_pts(ns+nm+1:end,:)=spotBG.fitRes;
            end
    end
else
    final_pts=[];
end
if verbose
    if nSingle-ns>0
        disp(['Removed Singles that do not fit: ',num2str(nSingle-ns)]);
    end
    if nMultiple-nm>0
        disp(['Removed multiples that are too dim: ',num2str(nMultiple-nm)]);
    end
    if nBG-nb>0
        disp(['Removed bg spots that are too dim: ',num2str(nBG-nb)]);
    end
    disp(['uLocalizeFitCC, Time spent: ', num2str(toc)]);
end
detection_res=struct('single',spotSingle,'multiple',spotMultiple,'spotBG', spotBG,'cc',cc);
end

