function p=uLocalizeInitPara(varargin)
%Obj: initialize the parameter structure for the uLocalize
%To change the default parameters, use the ParameterName, Value pair
%For example:
%   p=initULocalizePara()
%   p=initULocalizePara('sigma_xy',1.1, 'sigma_z',1.2, 'thickness',2,'method','GaussianMask');
ip=inputParser;
ip.addParameter('numdim',3,@(x) ismember(x,[2,3]));
ip.addParameter('sigma_xy',1, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('sigma_z',1, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('cutsize',3, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('TS_cutsize', 4, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('thickness',1, @(x) isnumeric(x) && isscalar(x) && x>0);

ip.addParameter('threshUnit', 'sd', @(x) ismember(lower(x), {'sd','absolute','automatic'}));
ip.addParameter('threshLevel',3, @(x) isnumeric(x) && isscalar(x) && x>0);

ip.addParameter('tol',0.001,@(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('max_count',50000,@(x) isnumeric(x) && isinteger(x) && x>0);
ip.addParameter('ROISizeOption','small',@(x) ismember(upper(x),{'small','large'})); %whether fittng window contains the thickness part or not, small: not; large: yes
ip.addParameter('bgCorrMethod', 'plane', @(x) ismember(lower(x),{'plane','quad', 'median', 'local plane', 'local quad', 'local median'})); %the local is the same as the other one
ip.addParameter('NormByCellLimit', 10, @(x) isnumeric(x) && isscalar(x) && x>0); %The number of spot above which to consider normalization by cell
ip.addParameter('ISingleLowBound', 0.32, @(x) isnumeric(x) && isscalar(x) && x>0); %The normalized intensity below which it is rejected as background spots
ip.addParameter('InternalParallel', false, @(x) islogical(x) && isscalar(x)); % use parallel inside fitting one image

%filter parameters
ip.addParameter('filterMethod','LoGRaj',@(x) ismember(upper(x),{'BANDPASS','LOGRAJ','LOGFFT'}));
ip.addParameter('filterSigma',1.2,@(x) isnumeric(x) && x>0);  %the sigma value of the LoG fiter
ip.addParameter('filterWindowSize',5, @(x) isnumeric(x) && x>0); %the window size of the Log fiter
ip.addParameter('nlo',7, @(x) isnumeric(x) && isscalar(x) && x>0);  %the bandpass filter, low freq
ip.addParameter('nhi',3, @(x) isnumeric(x) && isscalar(x) && x>0);  %the bandpass filter, high freq

ip.addParameter('ROIsize',2, @(x) isnumeric(x) && isscalar(x) && x>0); %local maximum window size
ip.addParameter('method','GaussianMask', @(x) ismember(upper(x),{'GAUSSIANMASK','GAUSSIANFIT'}));
%Dense mode parameters
ip.addParameter('detectionMode','LM', @(x) ismember(upper(x),{'LM','CC','GRANULE'}));
%ip.addParameter('denseMode',false, @(x) islogical(x) && isscalar(x)); %I have removed this parameter and changed it to detectionMode
ip.addParameter('aMin',5, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('aMax',500, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('eMin',0, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('eMax',0.55, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('numDilation',3, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('ISum2IFitRatio',1, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('bg_extension',1, @(x) isnumeric(x));

ip.parse(varargin{:});
    p.numdim=ip.Results.numdim;
    p.sigma_xy=ip.Results.sigma_xy;
    p.sigma_z=ip.Results.sigma_z;
    p.cutsize=ip.Results.cutsize;    %
    p.cutwidth = [ceil(p.cutsize*p.sigma_xy), ceil(p.cutsize*p.sigma_z)]; %The actual cropped image size 2*cutwidth+1;
    if p.cutwidth(1)<2, p.cutwidth(1)=2; end
    if p.cutwidth(2)<2, p.cutwidth(2)=2; end
    p.TS_cutsize=ip.Results.TS_cutsize;
    p.TS_cutwidth=[ceil(p.TS_cutsize*p.sigma_xy), ceil(p.TS_cutsize*p.sigma_z)];
    if p.TS_cutwidth(1)<3, p.cutwidth(1)=3; end %set TS_cutwidth between [3, 7]; prefer 5
    if p.TS_cutwidth(2)<3, p.cutwidth(2)=3; end
    if p.TS_cutwidth(1)>7, p.cutwidth(1)=7; end
    if p.TS_cutwidth(2)>7, p.cutwidth(2)=7; end
    p.thickness=ip.Results.thickness;      %With the new implementation, thickness could be set to 2
    
    p.thresh(1).level=ip.Results.threshLevel;
    p.thresh.unit = ip.Results.threshUnit;
    p.thresh.threshInt = p.thresh.level;
    
    p.tol=ip.Results.tol;
    p.max_count=ip.Results.max_count; %maximal number of points allowed
    p.bgCorrMethod = ip.Results.bgCorrMethod; % the local background correction method: can be plane, median, quad, see fitBGImgWrapper
    p.NormByCellLimit = ip.Results.NormByCellLimit; %number of points beyond which to normalize by cell, this actually should be with FQ version only
    p.ISingleLowBound = ip.Results.ISingleLowBound; %The normalized intensity below which it is rejected as background spots
    p.InternalParallel = ip.Results.InternalParallel; %use parallel in fitting singles
    
    p.filter.filterMethod = ip.Results.filterMethod;
    p.filter.nlo=ip.Results.nlo;    %The low frequency component: default 10
    p.filter.nhi=ip.Results.nhi;    %The high freqency component: default 3
    p.filter.numdim=2;  
    p.filter.width=0.01;
    p.filter.filterSigma=ip.Results.filterSigma;
    p.filter.filterWindowSize=ip.Results.filterWindowSize;
    if p.filter.filterWindowSize<5
        p.filter.filterWindowSize=max(2*round(2*p.filter.filterSigma)+1,5);
    end
    p.maxIter=100;  %for the Gaussian Mask iteration
    p.ROIsize=ip.Results.ROIsize;    %The distance between neiboring local maximum used in the pre-detection
    p.method=ip.Results.method; %Possible: 'GaussianFit3D', 'GaussianFit2D', 'GaussianMask2D'
    %The dense mode parameters
    p.detectionMode=ip.Results.detectionMode;
    p.denseMode=false;
    p.ROISizeOption = ip.Results.ROISizeOption; %whether fittng window contains the thickness part or not, small: not; large: yes
    p.aMin=ip.Results.aMin;
    p.aMax=ip.Results.aMax;
    p.eMin=ip.Results.eMin;
    p.eMax=ip.Results.eMax;
    p.numDilation=ip.Results.numDilation;
    p.ISum2IFitRatio=ip.Results.ISum2IFitRatio;
    p.bg_extension=ip.Results.bg_extension;
    
    p.dx=108.3;
    p.dy=108.3;
    p.dz=300;
end