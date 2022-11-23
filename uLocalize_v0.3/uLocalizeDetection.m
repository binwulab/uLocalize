function [final_pts, p, smooth, detection_res, fgImg]=uLocalizeDetection(img, p, maskType, mask, varargin)
%Objective: Perform spot detection using Gaussian Mask or Gaussian Fitting
%algorithm. 
%Input: 
%       img: a single 2D or 3D image
%       p: the parameter
%       maskType: 0: use the std of the whole image; 
%             1: each cell uses its own std; Default
%             2: all cells use the same std
%       mask: mask to define the region of interest
%       varargin: the (ParameterName,ParameterValue) pairs, including 
%               'method',default,p.method, will be overriden
%               'thresh', default, p.thresh
%               'verbose', default,1
%               'filteredImage', filtered Image, to save time
%               'forceFilter', logical, to force the program to filter
%Output:
%   final_pts: 
%       for GaussianMask2D: [yc,xc,N0,Chi2]
%       for GaussianFit2D: [yc,xc,N0,sigma_xy]
%       for GaussianMask3D: [yc,xc,zc,N0,err]
%       for GaussianFit3D: [yc,xc,zc,N0,sigma_xy,sigma_z]
%   p: the parameter structure, to pass back changed parameter here
%   smooth: to save the filtered smooth image because it is time consuming
%       to calculate it
%   detctionResults: depending on the detectionMode, the meaning is different
%       LM: detection_res=struct('rough_pts', rough_pts, 'final_pts', final_pts);
%       CC: detection_res=struct('single',spotSingle,'multiple',spotMultiple,'cc',cc);
%       GRANULE: detection_res=cc;
%   fgImg: The foreground Image, each pixel containing the number of object
%       that pixel belongs to. Useful for subsequent analysis such as TS 
%       
%   
%History   
%   B.W., Jan, 2021
%   B.W., Feb, 2021, change bgImg to fgImg
%   B.W., Mar, 2021, change the parameter denseMode to detectionMode,
%       adding the granule detection module

ip=inputParser;
ip.addRequired('img',@(x) isnumeric(x) && ~isscalar(x));
ip.addRequired('p',@(x) isstruct(x));
ip.addOptional('maskType',1,@(x) ismember(x,0:2));
ip.addOptional('mask',[],@(x) isempty(x) || ~isscalar(x));
ip.addParameter('method',p.method, @(x) ismember(x,{'GaussianMask','GaussianFit'}));
%ip.addParameter('thresh',p.thresh,@(x) isstruct(x));
ip.addParameter('filteredImage',[],@(x) isnumeric(x) && ~isscalar(x));
ip.addParameter('forceFilter', false, @(x) islogical(x) && isscalar(x));
ip.addParameter('calibration', false, @(x) islogical(x) && isscalar(x));
ip.addParameter('verbose',true, @(x) isnumeric(x) | islogical(x));

ip.KeepUnmatched = true; %So it is OK if there is unmatched Parameter-Value pair
ip.parse(img,p, maskType, mask, varargin{:});

sz=size(img);
p.numdim=numel(sz);
maskType=ip.Results.maskType;
if isempty(ip.Results.mask) || maskType==0 %If mask is undefined, set to the whole image
    mask=ones(sz(1:2)); %Only make the 2D mask
end
%thresh=ip.Results.thresh;
p.method=ip.Results.method;     %Override the method in p. But it will not pass back since p is not in the output argument
verbose=ip.Results.verbose;

%% smooth the image with band pass filter
if isempty(ip.Results.filteredImage) || ip.Results.forceFilter
    smooth=uLocalizeFilter(img,p);
else
    smooth=ip.Results.filteredImage;
end

switch upper(p.detectionMode) 
    case 'CC'   %
        %% Predetection with Connected component thresholding
        cc=predetectCC(smooth, p.thresh, maskType, mask, [p.aMin, p.aMax],[p.eMin,p.eMax], verbose);
        if nargout<3
            clear('smooth');    %release the memory,if not saved for later
        end

        %% Detection objects one by one, returned fgImg
        [final_pts, p, detection_res, fgImg]=uLocalizeFitCC(img,p,cc,ip.Results.calibration,verbose); 
    case 'LM'
        %% Predetection with local maximum algorithm
        rough_pts=predetectLM(smooth, p.thresh, maskType, mask, p.ROIsize, p.max_count, verbose);
        if nargout<3
            clear('smooth');    %release the memory,if not saved for later
        end

        %% Detection points one by one, returned fgImg
        [final_pts,p, detection_res, fgImg]=uLocalizeFitLM(img,p,rough_pts,verbose); 
    case 'GRANULE'  %granule detection
        cc=predetectGranule(smooth, p.thresh, maskType, mask, p.aMin, p.eMax, verbose);
        if nargout<3
            clear('smooth');
        end
        %Note the detection result is the cc structure
        [final_pts, p, detection_res,fgImg]=uLocalizeFitGranule(smooth, p, cc, verbose); %No te I am using the smoothed image for this. To really quantify, we need to use the original image
        %fgImgOrBW=dummy;    %This is unused at the moment, just to satisfy the calling semantics
    otherwise
        disp('uLocalizeDetection: unknown detection mode');
        return
end
end