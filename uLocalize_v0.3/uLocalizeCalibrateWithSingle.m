function p=uLocalizeCalibrateWithSingle(img, p, cc, bw)
%Objective: Calibrate parameters for calculating and interpreting
%integrated intensities, using well defined singlets
%Input: 
%   img: the raw image
%   p: the detection parameter structures
%   cc: the connected component using thresholded predetection
%   bw: the mask file containing all the thresholded spots
%Output
%   p: the output detection parameters with added/updated parameters
%       bg_extension: the average number of pixel from the obj to the boundary of the cropped image
%       numDilation: the times of dilation needed to get the intensity sum equal to the fitted intensity
%       ISum2IFitRatio: the ratio between ISum and IFit, used to convert ISum to IFit for complext objects
%History
%   BW, Jan 2021

%Fit all single candidates
numdim=ndims(img);
nSingle=cc.nObjects(1);
ObjCtr=cc.ObjCtr(1:nSingle,:);
pixelIdxList=cc.PixelIdxList(1:nSingle);   %This results in better performance for parallel
[fitRes, fitStatus, ISum, bg_mean, bg_sd, nPix, bg_ext] = ...
    uLocalizeFitSingleWrapper(img, bw, p, ceil(ObjCtr),pixelIdxList);

%[fitRes, fitStatus, ISum,bg_mean,bg_sd,nPix, bg_ext] = uLocalizeFitSingleWrapper(img, p, cc, bw);
indGood= find(fitStatus==1);    %choose those spots with good fit
if numel(indGood) < 1   %no good spots
    disp('uLocalizeCalibrateWithSingle: too bad, cannot find a good spot');
    return;
end
fitRes=fitRes(indGood,:);   %only use good spots
ISum=ISum(indGood);
IFit=fitRes(:,numdim+1);
ISingleMed=median(IFit);    %Median fit intensity
avgIFit=median(IFit);         %mean of fit intensity
avgISum=median(ISum);         %mean of intensity sum
if avgISum < avgIFit        %The thresholded spot is too small to calculate intensity, needs to grow
    disp(['Spot is too small: IFit=', num2str(median(IFit)), '; ISum=', num2str(median(ISum))]);
    growSpots=true;
    %pTmp=p;
else
    disp('uLocalizeCalibrateWithSingle: the spot size is fine, no need to dilate');
    growSpots=false;
end
numDilation=0;
while growSpots
    numDilation=numDilation+1;
    %Note: dilate from the original image, not from the previous dilated
    %one. And it dilate all objects
    [fgImg, ccTmp, areaDiff] = dilateObjects(bw, cc, numDilation);
    ObjCtr=ccTmp.ObjCtr(1:nSingle,:);
    pixelIdxList=ccTmp.PixelIdxList(1:nSingle);   %This results in better performance for parallel
    [fitRes, fitStatus, ISum, bg_mean, bg_sd, nPix, bg_ext] = ...
        uLocalizeFitSingleWrapper(img, fgImg, p, ceil(ObjCtr),pixelIdxList);
    %[fitRes, fitStatus, ISum,bg_mean,bg_sd,nPix, bg_ext] = uLocalizeFitCCSingle(img, p, ccTmp, fgImg);
    indGood= find(fitStatus==1);    %only select good spots
    fitRes=fitRes(indGood,:);
    IFit=fitRes(:,numdim+1);
    ISingleMed=median(IFit);
    ISum=ISum(indGood);
    bg_sd=bg_sd(indGood);
    avgISum2=median(ISum);
    dI=avgISum2-avgISum;    %Note: dI is the differential increment from successive dilation
    %Should use only the good single spots. 
    %Note: the objects are arranged in the order of [single, multiple, zeros]. So it is ok to do this
    areaDiff=median(areaDiff(indGood)); 
    %pTmp.amax=pTmp.amax+areaDiff;
    if dI <=2*median(bg_sd)*sqrt(areaDiff) || avgISum2 > ISingleMed %areaDiff is the changing number of pixel
        growSpots = false;
    else
        avgISum=avgISum2;
    end
    disp(['Iteration ',num2str(numDilation), '; IFit=', num2str(median(IFit)), ...
        '; ISum=', num2str(median(ISum)), '; dI=', num2str(dI), '; noise=', num2str(2*median(bg_sd)*sqrt(areaDiff))]);
end
p.bg_extension=mean(bg_ext,1); %average the bg_ext for different spot
p.ISum2IFitRatio=median(ISum)/ISingleMed;  %The calibration ratio between intensity sum and fitted intensity
p.numDilation=numDilation;
if numDilation>1
    p.aMax=p.aMax+areaDiff;
end
end