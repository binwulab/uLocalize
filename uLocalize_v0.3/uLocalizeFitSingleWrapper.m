function [fitRes, fitStatus, ISum, bg_mean, bg_sd, nPix, bg_ext]= ...
    uLocalizeFitSingleWrapper(img, fgImg, p, ObjCtrs, pixelIdxList)
%Objective: this function was originally in the uLocalizeFitCCSingle, taken
%out for parallization and quantification of bg spots. Do the fitting for a
%single spot, and generalize it to the LM mode. So we have a unified single
%spot fitter
%Input:
%   image: the raw image
%   fgImg: the foreground image of the objects. For final quantification, each pixel
%       represents the number of object containing it
%   p: the detection parameter structure
%   ObjCtr: the center of the object, or the local maximum pixel
%   pixelIdxList: the pixel list of the object in the cc structure
    %   cc: the connected component structure of candidate objects with extra fields
    %       Connectivity: 4 for 2D, 18 for 3D
    %       ImageSize: the raw image size
    %       NumObjects: the number of objects
    %       nObjects = [nSingle, nMultiple, nZeros], nSingle is the number of
    %           single candidates, nMultiple: the number of complex obj, nZeros:
    %           the number of background noisy spot
    %       PixelIdxList: array linear indices of each objects, it should be
    %           sorted in the order of [Single, Multiple, Zeros]
%History: 
%   B.W. Jan 2021
%   B.W. Aug 2021: unify the single spot fitter for both LM and CC mode

numdim=ndims(img);
nd=size(ObjCtrs,1); %number of points
if nargin<5
    pixelIdxList=cell(nd,1);   %Note this is required to call the nested function
end

if strcmpi(p.method,'GaussianFit')
    fitRes=zeros(nd, numdim*2);    %[y,x, z,I,sigma_xy, sigma_z] for 3D [y,x,I,sigma_xy]
else
    fitRes=zeros(nd, numdim+2);    %[x,y,z,I,err0]
end
if nargout>1
    fitStatus=zeros(nd,1);  %0:N/A; 1:Good; 2:diverge; 3:center not in the obj
    ISum = zeros(nd,1);     %The sum of intensity under the object
    bg_mean = zeros(nd,1);  %The background mean value of the interpolated plane
    bg_sd =zeros(nd,1);     %the std of the surrounding pixel
    nPix = zeros(nd,1);     %num of pixel for the object
    bg_ext = zeros(nd,numdim-1);    %2D: [ext_xy], 3D: [ext_xy, ext_z]
end
if p.InternalParallel
    parfor i=1:nd
        %feval(funcHandle, i);
        [fitRes(i,:), fitStatus(i), ISum(i), bg_mean(i), bg_sd(i), nPix(i), bg_ext(i,:)] = ...
            fitSingleFunc(img, fgImg, p, ObjCtrs(i,:), pixelIdxList{i});
    end
else
    for i=1:nd
%         i
%         if i==3002
%             dummy=input('hit me');
%         end
        %feval(funcHandle, i);
        [fitRes(i,:), fitStatus(i), ISum(i), bg_mean(i), bg_sd(i), nPix(i), bg_ext(i,:)] = ...
            fitSingleFunc(img, fgImg, p, ObjCtrs(i,:), pixelIdxList{i});
    end
end
%% Implementation with nested function. Somehow did not work with parallel
% funcHandle = @fitSingleFuncNested;
% if p.InternalParallel
%     parfor i=1:nd
%         feval(funcHandle, i);
%     end
% else
%     for i=1:nd
%         feval(funcHandle, i);
%     end
% end


    function fitSingleFuncNested(indx)        
        [imgCorr,fgImg2,ObjCtr2, bndBox, ~, bgPlane_mean, bgCorr_sd, pixelIdxList2, bgExt]= ...
            genBGCorrForPointByLinInterp(img, fgImg, ObjCtrs(indx,:),  p.cutwidth, p.thickness, p.ROISizeOption, p.bgCorrMethod, pixelIdxList{indx});
        
        if strcmpi(p.detectionMode,'CC') && ~isempty(pixelIdxList{indx})   %This would be in CC mode
            currObjPixels=fgImg2(pixelIdxList2);   %Note fgImg2 is not binary, so needs to store the value. After fitting the current spot, put it back
            %Note we have set the current obj to 0 in the cropped mask bw2. This is
            %important because in the cropped image, there may be still pixels
            %contained in other objects. These pixels are discarded during fitting
            fgImg2(pixelIdxList2)=0;
        else
            fgImg2=fgImg2-1;    %This would be LM mode
        end
        
        switch upper(p.method)
            case 'GAUSSIANFIT'  %To be implemented
                [pt, status]=gaussianFitWrapper(imgCorr, ObjCtr2, 'integrated gaussian', fgImg2);
                if status==0
                    fitRes(indx,:)=pt(1:numdim*2);
                else
                    fitRes(indx,:)=nan;
                end
            otherwise   %The default is GaussianMask
                [pt, status]=gaussianMaskWrapper(imgCorr, ObjCtr2, [p.sigma_xy, p.sigma_z], p.maxIter, p.tol,fgImg2);
                if status==0
                    fitRes(indx,:)=pt(1:numdim+2);
                else
                    fitRes(indx,:)=nan;
                end
        end
        fitRes(indx,1:numdim)=pt(1:numdim)+bndBox(1,1:numdim)-0.5; %change it to the original coordinate, -1.5 is usd to offset the indexing
        
        if strcmpi(p.detectionMode,'CC') && ~isempty(fitStatus)   %CC mode and output is requested
            fgImg2(pixelIdxList2)=currObjPixels;    % restore the value
            %fit status:
            %0: n/a
            %1: OK, Note the default was set to 1
            %2: gaussian mask diverged
            %3: center of the spot not within the Object
            ctr=pt(1:numdim);
            ind=sub2ind2(size(imgCorr),ceil(ctr)); %convert the obj center to linear index in the cropped image
            fitStatus(indx)=1;
            if ctr(1) == -1     %not converged
                fitStatus(indx)=2;
            elseif numel(find(ind==pixelIdxList2))==0  %If the object center is not within the object
                fitStatus(indx)=3;
            end
            
            %calculate the intensity sum, Note in bw2, each pixel value represents how many obj contained that pixel
            %INormByObjNum = imgCorr(pixelIdxList2)./double(bw2(pixelIdxList2));
            ISum(indx)=sum(imgCorr(pixelIdxList2)./double(fgImg2(pixelIdxList2)));
            bg_mean(indx)=bgPlane_mean;
            bg_sd(indx)=bgCorr_sd;
            nPix(indx)=numel(pixelIdxList2);
            bg_ext(indx,:)=bgExt;
        end
    end
end


function [fitRes,fitStatus,ISum, bg_mean, bg_sd, nPix, bg_ext] = fitSingleFunc(img, fgImg, p, ObjCtr, pixelIdxList)
numdim=ndims(img);
if strcmpi(p.method,'GaussianFit')
    fitRes=zeros(1, numdim*2);    %[y,x, z,I,sigma_xy, sigma_z] for 3D [y,x,I,sigma_xy]
else
    fitRes=zeros(1, numdim+2);    %[x,y,z,I,err0]
end

[imgCorr,fgImg2,ObjCtr2, bndBox, ~, bgPlane_mean, bgCorr_sd, pixelIdxList2, bgExt]= ...
    genBGCorrForPointByLinInterp(img, fgImg, ObjCtr,  p.cutwidth, p.thickness, p.ROISizeOption, p.bgCorrMethod, pixelIdxList);

if strcmpi(p.detectionMode,'CC') && ~isempty(pixelIdxList)   %This would be in CC mode
    currObjPixels=fgImg2(pixelIdxList2);   %Note fgImg2 is not binary, so needs to store the value. After fitting the current spot, put it back
    %Note we have set the current obj to 0 in the cropped mask bw2. This is
    %important because in the cropped image, there may be still pixels
    %contained in other objects. These pixels are discarded during fitting
    fgImg2(pixelIdxList2)=0;
else
    fgImg2=fgImg2-1;    %This would be LM mode
end

switch upper(p.method)
    case 'GAUSSIANFIT'  %To be implemented
        [pt, status]=gaussianFitWrapper(imgCorr, ObjCtr2, 'integrated gaussian', fgImg2);
        if status==0
            fitRes(1,:)=pt(1:numdim*2);
        else
            fitRes(1,:)=nan;
        end
    otherwise   %The default is GaussianMask
        [pt, status]=gaussianMaskWrapper(imgCorr, ObjCtr2, [p.sigma_xy, p.sigma_z], p.maxIter, p.tol,fgImg2);
        if status==0
            fitRes(1,:)=pt(1:numdim+2);
        else
            fitRes(1,:)=nan;
        end
end
fitRes(1,1:numdim)=fitRes(1:numdim)+bndBox(1,1:numdim)-0.5; %change it to the original coordinate, -1.5 is usd to offset the indexing

if strcmpi(p.detectionMode,'CC') && ~isempty(pixelIdxList)    %CC mode and output is requested
    fgImg2(pixelIdxList2)=currObjPixels;    % restore the value
    %fit status:
    %0: n/a
    %1: OK, Note the default was set to 1
    %2: gaussian mask diverged
    %3: center of the spot not within the Object
    ctr=pt(1:numdim);
    ind=sub2ind2(size(imgCorr),ceil(ctr)); %convert the obj center to linear index in the cropped image
    fitStatus=1;
    if ctr(1) == -1     %not converged
        fitStatus=2;
    elseif numel(find(ind==pixelIdxList2))==0  %If the object center is not within the object
        fitStatus=3;
    end
    
    %calculate the intensity sum, Note in bw2, each pixel value represents how many obj contained that pixel
    %INormByObjNum = imgCorr(pixelIdxList2)./double(bw2(pixelIdxList2));
    ISum=sum(imgCorr(pixelIdxList2)./double(fgImg2(pixelIdxList2)));
    bg_mean=bgPlane_mean;
    bg_sd=bgCorr_sd;
    nPix=numel(pixelIdxList2);
    bg_ext=bgExt;
else %still need to return to keep the calling syntacs. 
    fitStatus=nan;
    ISum=nan;
    bg_mean=nan;
    bg_sd=nan;
    nPix=nan;
    bg_ext=nan;
end
end