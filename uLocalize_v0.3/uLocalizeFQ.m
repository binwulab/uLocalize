function [cell_pts, cell_ts, cell_prop, final_pts,smooth,p, cell_granule]=uLocalizeFQ(img, p, outlineFile, outputFile, maskType,  cellRegion,  varargin)
%%Objective:
% Call the uLocalize with the outline defined by the FQ. 
%Input:
%       img: the 3D stack of FISH image
%       p: the parameter structure for uLocalize
%       outlineFile: the full path of the outline file
%       Optional input parameters
%           outputFile: the name of the output .txt File, read by FQ
%           maskType: 
%               0: Full image
%               1: calculate threshold for each cell individually
%               2: calculate the threshold for all cells in the image
%           cellRegion: 
%               'cytosol': cytoplasmic only
%               'nuclear': nuclear only
%               'cell': the whole cell
%       Optional Parameter-value pair
%           'method': 'GaussianMask', 'GaussianFit'
%           'threshold': threshold used for predetection in units STD
%           'filteredImage': smoothed 3D stack 
%           'forceFilter':ignore the provided smoothed image and filter
%           'calibration': used for the dense CC detection mode
%           'granuleCC': the granule CC structure 
%Output
%       cell_pts: a cell array containing the final points for each cell
%       cell_ts: a cell array containg the final TS for each cell
%                %[y,x,z,I,sxy,sz,bg,INormByCell],     
%       cell_prop: the cell property structure for FQ
%       final_pts: the final_pts containing all spots
%               if GaussianFit3D: [yc, xc, zc, Int, sigma_xy, sigma_z, bg]
%               if GaussianMask3D: [yc, xc, zc, Int, err];
%       smooth: the smoothed image, which can be saved and loaded again to save future run time
%       p: the parameter structure, to pass back the changed parameter here
%       cell_granule: a cell array containing the granule structure for each cell
%Update History   
%   B.W., Jan, 2021
%   B.W., Mar, 2021, added the granule detection mode and granule quantification
%   B.W., Aug, 2021, 
        %update the CC mode to include BG spots. 
        %try to unify the pixel offset from single spot, ts, and granule

ip=inputParser;
ip.addRequired('img', @(x) isnumeric(x) && ~isscalar(x));
ip.addRequired('p', @(x) isstruct(x));
ip.addRequired('outlineFile',@(x) isstr(x) && exist(x,'file')==2);
ip.addOptional('outputFile',[],@(x) isstr(x));
ip.addOptional('maskType',1, @(x) ismember(x,0:2));
ip.addOptional('cellRegion','cytosol',@(x) ismember(x,{'cytosol','nuclear','cell'}));
ip.addParameter('method',p.method, @(x) ismember(x,{'GaussianMask','GaussianFit'}));
% ip.addParameter('thresh',p.thresh,@(x) isstruct(x));
ip.addParameter('filteredImage',[],@(x) isnumeric(x) && ~isscalar(x));
ip.addParameter('forceFilter', false, @(x) islogical(x) && isscalar(x));
ip.addParameter('calibration', false, @(x) islogical(x) && isscalar(x));    %For dense mode only
ip.addParameter('granuleCC', [], @(x) isstruct(x) && isfield(x,'PixelIdxList')); %input the granule objects
ip.addParameter('verbose',false, @(x) isnumeric(x) | islogical(x));

ip.parse(img,p,outlineFile,outputFile,maskType,cellRegion,varargin{:});

outputFile=ip.Results.outputFile;
maskType=ip.Results.maskType;
cellRegion=ip.Results.cellRegion;
p.method=ip.Results.method; %Override the method in p. But it will not pass back since p is not in the output argument
% p.thresh=ip.Results.thresh;
numdim=ndims(img);

if strcmpi(p.detectionMode, 'GRANULE') && ~isempty(ip.Results.granuleCC)
    disp('uLocalizeFQ: the granule detection mode and granuleCC optional parameter should not be given at the same time');
    return;
end

%Read in the cell outlines and define the region of quantification
[cell_prop,par_microscope,file_names,flag_file,version,size_img,comment]=FQ_load_results_WRAPPER_v2(outlineFile,'');
numCells=numel(cell_prop);
sz=size(img);
ny=sz(1);
nx=sz(2);
mask=zeros(ny,nx,'uint16');
for iC=1:numCells
    bw=genCellRegion(cell_prop(iC), sz, cellRegion) & ~mask;    %generate the ROI for the current cell, depending on the cellRegion
    if (~isempty(cell_prop(iC).pos_TS)) %Remove TS for normal quantification
        bwTS=false(ny,nx);
        pos_TS=cell_prop(iC).pos_TS;
        for j=1:numel(pos_TS)
            bwTS=bwTS | poly2mask(pos_TS(j).x, pos_TS(j).y, ny,nx); %remove each TS
        end
        bw = bw & (~bwTS);
    end
    mask=mask+uint16(bw)*iC; %label the region with the cell number
end

%If granuleCC is given, remove the these ROIs from the subsequent FISH spot detection
if ~isempty(ip.Results.granuleCC)
    granuleCC=ip.Results.granuleCC; 
    for iC=1:numel(granuleCC) %If granuleCC is an array, it means it is for different cells
        for i=1:granuleCC(iC).NumObjects
            %Note if we treat these granule as TS, it is equivalent to
            %remove the bounding box. If not, these granules will be removed
            %here. 
            %remove the whole bounding box
            bndBox=getObjBndBox(granuleCC(iC).PixelIdxList{i},granuleCC(iC).ImageSize, p.bg_extension);
            mask(bndBox(1,1):bndBox(2,1), bndBox(1,2):bndBox(2,2))=0;
            %or remove granule pixel + one dilation layer
            %mask(dilateOneObject(granuleCC(iC).PixelIdxList{i}, granuleCC(iC).ImageSize,1))=0;
        end
    end
end

%calling the uLocalize
[final_pts, p,smooth,detection_res,fgImgOrBW]=uLocalizeDetection(img,p,maskType,mask,varargin{:});

%Convert the results to FQ format
if ~isempty(final_pts)
    [cell_prop,cell_pts]=convertULoc2FQ(final_pts, sz, p, cellRegion, cell_prop, par_microscope);
    parameters=struct('path_save','', 'cell_prop', cell_prop, ...
        'par_microscope', par_microscope, 'file_names', file_names,...
        'version', version, 'flag_type', 'spots');
    if ~isempty(outputFile)
        FQ_save_results_v1(outputFile, parameters);
    end
else
    cell_pts=[];
end

%Now Handle the transcription sites
cell_ts=cell(1,numCells);
IMedAll=median(final_pts(:,numdim+1));  %median intensity for spots in the whole image
for iC=1:numCells
    if (~isempty(cell_prop(iC).pos_TS)) %Remove TS for normal quantification
        pos_TS=cell_prop(iC).pos_TS;
        ts_tmp=zeros(numel(pos_TS),2*numdim+1);  %[y,x,z,I,sxy,sz,bg,INorm]
        IMedByCell=IMedAll;     %If there are too little spots in the cell, use the median of the whole image
        if size(cell_pts{iC},1)>  p.NormByCellLimit    %only if there are more than 10 spots in the cell, we define IMedByCell
            IMedByCell = median(cell_pts{iC}(:,numdim+1));
        end
        for j=1:numel(pos_TS)
            yc=ceil(mean(pos_TS(j).y));
            xc=ceil(mean(pos_TS(j).x));
            cutwidth_xy=ceil((max(pos_TS(j).x)-min(pos_TS(j).x))/2);
            cutwidth_z=ceil(p.cutwidth(2)*cutwidth_xy/p.cutwidth(1));    %scale the cutwidth_z accordingly
            cutwidth=[cutwidth_xy cutwidth_z];
            %cutwidth=p.TS_cutwidth;
            thickness=p.thickness;
            [int,indx]=max(img(:));
            [~,~,zc]=ind2sub(size(img),indx);   %find the maximal z-plane
            [int,zc]=max(img(yc,xc,:));
            pt0=double([yc, xc, zc, int]);
            [imgCorr, fgImg2, pt2, bndBox]=genBGCorrForPointByLinInterp(img, fgImgOrBW, pt0, cutwidth, thickness, 'small','quad'); %generate a cropped, local background corrected image
            %Note: if the TS is masked out in the beginning, fgImg2 should
            %be 0 around the TS. We've set fgImg2 to be uint16. So within
            %the box, it should be 0. 
            pt=gaussianFitWrapper(imgCorr, pt2, 'integrated gaussian',fgImg2);
            ts_tmp(j,1:2*numdim)=pt(1:2*numdim);
            ts_tmp(j,1:numdim)=ts_tmp(j,1:numdim)+bndBox(1,1:numdim)-0.5;  % Add the offset to the position 2021_08_15
            ts_tmp(j,2*numdim+1) = ts_tmp(j,numdim+1)/IMedByCell;       %normalize by intensity in the current cell, 2021_08_15, Note this number is over-written in the CC mode by only considering singles
        end
        cell_ts{iC}=ts_tmp;
    end
end



%% If the detection mode is granule, need to distribute the granule into cells. Note the detection_res is supposed to be the cc structure
if strcmpi(p.detectionMode, 'GRANULE')   
    granuleCtr=detection_res.ObjCtr;    %detection_res is cc here, should be one element
    ind=ny*(ceil(granuleCtr(:,2))-1)+round(granuleCtr(:,1));  %converts granule center to linear index
    cell_granule=repmat(detection_res, numCells,1);
    for iC=1:numCells
        bw=genCellRegion(cell_prop(iC), sz, cellRegion); %generate the mask for the current cell
        indSel=find(bw(ind));   %Note: bw(ind) will be the same length as fina_pts, with 0 or 1. So this returns the indices of points in the current cell
        cell_granule(iC).PixelIdxList=detection_res.PixelIdxList(indSel);
        cell_granule(iC).NumObjects=numel(cell_granule(iC).PixelIdxList);
        cell_granule(iC).ObjCtr=detection_res.ObjCtr(indSel,:);
        cell_granule(iC).ISum=detection_res.ISum(indSel);
    end
end

%% If the detection mode is CC, need to distribute the detection results into cells too. 
%   This was the definition in uLocalizeCC: detection_res=struct('single',spotSingle,'multiple',spotMultiple,'cc',cc);
if strcmpi(p.detectionMode, 'CC')
    spotSingle=detection_res.single;
    spotMultiple=detection_res.multiple;
    spotBG=detection_res.spotBG;
    
    singleCtr=spotSingle.fitRes;
    if ~isempty(singleCtr)
        indSingle=ny*(ceil(singleCtr(:,2))-1)+round(singleCtr(:,1));  %converts single center to linear index
    else
        indSingle=[];
    end
    multipleCtr=spotMultiple.fitRes;
    if ~isempty(multipleCtr)
        indMultiple=ny*(ceil(multipleCtr(:,2))-1)+round(multipleCtr(:,1));  %converts multiple center to linear index
    else
        indMultiple=[];
    end
    bgCtr=spotBG.fitRes;
    if ~isemtpy(bgCtr)
        indBG=ny*(ceil(bgCtr(:,2))-1)+round(bgCtr(:,1));  %converts multiple center to linear index
    else
        indBG=[];
    end
    
    ISingleMedAll=median(spotSingle.fitRes(:,numdim+1));  %median intensity for singles in the whole image
    ISingleQuartileAll=quantile(spotSingle.fitRes(:,numdim+1), 0.25);
    for iC=1:numCells
        bw=genCellRegion(cell_prop(iC), sz, cellRegion); %generate the mask for the current cell
        indSel=find(bw(indSingle));   %Note: bw(ind) will be the same length as spotSingle, with 0 or 1. So this returns the indices of points in the current cell
        ISingleMedByCell=ISingleMedAll;     %If there are too little spots in the cell, use the median of the singles in the whole image
        ISingleQuartileByCell=ISingleQuartileAll;
        if numel(indSel)>  p.NormByCellLimit    %only if there are more than 10 spots in the cell, we define IMedByCell
            ISingleMedByCell = median(spotSingle.fitRes(indSel,numdim+1));
            ISingleQuartileByCell=quantile(spotSingle.fitRes(indSel,numdim+1), 0.25);
        end
        if numel(indSel)>0
            cell_prop(iC).single=struct('fitRes',spotSingle.fitRes(indSel,:), ...
                'INorm',spotSingle.INorm(indSel,:),     ...
                'fitStatus', spotSingle.fitStatus(indSel,:),    ...
                'ISum',spotSingle.ISum(indSel,:),...
                'bg_mean',spotSingle.bg_mean(indSel,:),     ...
                'bg_sd',spotSingle.bg_sd(indSel,:),     ...
                'nPix',spotSingle.nPix(indSel,:),   ...
                'bg_ext',spotSingle.bg_ext(indSel,:));
        else
            cell_prop(iC).single=struct('fitRes',[], ...
                'INorm',[],     ...
                'fitStatus', [],    ...
                'ISum',[],...
                'bg_mean',[],     ...
                'bg_sd',[],     ...
                'nPix',[],   ...
                'bg_ext',[]);
        end
        indSel=find(bw(indMultiple));   %Note: bw(ind) will be the same length as spotMultiple, with 0 or 1. So this returns the indices of points in the current cell
        if numel(indSel)>0
            INormByCell = spotMultiple.fitRes(indSel,numdim+1) /ISingleMedByCell;
            %nRNAByCell = max([(INormByCell> 1-ISingleQuartileByCell / ISingleMedByCell & INormByCell<1.5),round(INormByCell)],[], 2);
            nRNAByCell = max([(INormByCell> p.ISingleLowBound & INormByCell<1.5),round(INormByCell)],[], 2);
            cell_prop(iC).multiple=struct('fitRes',spotMultiple.fitRes(indSel,:),   ...
                'INorm', spotMultiple.INorm(indSel,:),  ...
                'ISum', spotMultiple.ISum(indSel,:),    ...
                'bg_mean',spotMultiple.bg_mean(indSel,:),   ...
                'bg_sd',spotMultiple.bg_sd(indSel,:), ...
                'nPix', spotMultiple.nPix(indSel,:),    ...
                'nRNA',spotMultiple.nRNA(indSel,:), ...
                'INormByCell', INormByCell,      ...
                'nRNAByCell',nRNAByCell);
        else
            cell_prop(iC).multiple=struct('fitRes',[],   ...
                'INorm', [],  ...
                'ISum', [],    ...
                'bg_mean',[],   ...
                'bg_sd',[], ...
                'nPix', [],    ...
                'nRNA',[], ...
                'INormByCell', [],      ...
                'nRNAByCell',[]);
        end
        indSel=find(bw(indBG));   %Note: bw(ind) will be the same length as spotSingle, with 0 or 1. So this returns the indices of points in the current cell
        if numel(indSel)>0
            INormByCell = spotBG.fitRes(indSel,numdim+1) /ISingleMedByCell;
            %nRNAByCell = max([(INormByCell>1-ISingleQuartileByCell/ISingleMedByCell & INormByCell<1.5),round(INormByCell)],[], 2);
            nRNAByCell = max([(INormByCell>p.ISingleLowBound & INormByCell<1.5),round(INormByCell)],[], 2);
            cell_prop(iC).spotBG=struct('fitRes',spotBG.fitRes(indSel,:), ...
                'INorm',spotBG.INorm(indSel,:),     ...
                'fitStatus', spotBG.fitStatus(indSel,:),    ...
                'ISum',spotBG.ISum(indSel,:),...
                'bg_mean',spotBG.bg_mean(indSel,:),     ...
                'bg_sd',spotBG.bg_sd(indSel,:),     ...
                'nPix',spotBG.nPix(indSel,:),   ...
                'bg_ext',spotBG.bg_ext(indSel,:), ...
                'nRNA', spotBG.nRNA(indSel,:),  ...
                'INormByCell', INormByCell,      ...
                'nRNAByCell',nRNAByCell);
        else
            cell_prop(iC).spotBG=struct('fitRes',[], ...
                'INorm',[],     ...
                'fitStatus', [],    ...
                'ISum',[],...
                'bg_mean',[],     ...
                'bg_sd',[],     ...
                'nPix',[],   ...
                'bg_ext',[], ...
                'nRNA', [],  ...
                'INormByCell', [],      ...
                'nRNAByCell',[]);
        end
        if ~isempty(cell_ts{iC})
            cell_ts{iC}(:,2*numdim+1) = cell_ts{iC}(:,numdim+1)/ISingleMedByCell;   %Note I re-calculate the normalized ts for the CC mode: 08-15-2021
        end
     end
end

%If granuleCC is given, it means the granule has been calcualted previously
%Now calculate the integrated intensity granule for this channel
if ~isempty(ip.Results.granuleCC)
    granuleCC=ip.Results.granuleCC;
    cell_granule=granuleCC;
    for i=1:numel(granuleCC)
        [ISum, dist2Ctr, dist2Obj]=calcIntSumGranule(img,p,granuleCC(i), cell_ts{i}); %this calculate the ISum for all granules in this cell
        cell_granule(i).ISum=ISum;    %Note: the ObjCtr should be identical, so only need to calculate the intensity
        cell_granule(i).dist2Ctr=dist2Ctr;
        cell_granule(i).dist2Obj=dist2Obj;
    end
end

end


function bw=genCellRegion(cell_prop, imageSize, cellRegion)
%For each cell given in the cell_prop, generate a mask with the same size
%of the original image. The region of the cell is determined by the string
%cell region: cytosol, nuclear, or the whole cell. 
imageSizeY=imageSize(1);
imageSizeX=imageSize(2);
roiX=cell_prop.x;
roiY=cell_prop.y;
bw=poly2mask(roiX,roiY,imageSizeY,imageSizeX);  %generate the ROI for the whole cell
if ~(isempty(cell_prop.pos_Nuc))    %only when nucleus is defined, can use the cellRegion
    nucROIX=cell_prop.pos_Nuc.x;
    nucROIY=cell_prop.pos_Nuc.y;
    bwNuc=poly2mask(nucROIX, nucROIY, imageSizeY, imageSizeX);  %generate the nuclear ROI

    switch lower(cellRegion)    %choose the cellular region
        case 'cytosol'  % cytoplasmic only
            bw = bw & (~bwNuc);
        case 'nuclear'  % nuclear
            bw=bwNuc;
        otherwise % 'cell'  %in this case, do nothing, just keep the whole cell
    end
elseif strcmpi(cellRegion, 'cytosol') ||  strcmpi(cellRegion, 'nuclear')
    disp('WARNING! uLocalizeFQ->genCellRegion: nucleus is not defined, cannot select cytosol or nuclear region, so choose the whole cell');
end
end


function [cell_prop,cell_pts]=convertULoc2FQ(final_pts, imageSize, p, cellRegion, cell_prop, par_microscope)
% convert the results obtained from uLocalize to the format for FQ
IMGSizeY=imageSize(1);
IMGSizeX=imageSize(2);
IMGSizeZ=imageSize(3);
pixSizeXY=par_microscope.pixel_size.xy;
pixSizeZ=par_microscope.pixel_size.z;
boundSizeXY=min([p.cutwidth(1),4]); %I don't know why, but FQ cannot handle more than 4, so I set this condition
boundSizeZ=min([p.cutwidth(2),4]);   
sigma_XY=p.sigma_xy;
sigma_Z=p.sigma_z;
numCells=numel(cell_prop);
cell_pts=cell(1,numCells);
ind=IMGSizeY*(ceil(final_pts(:,2))-1)+round(final_pts(:,1));  %converts all points to linear index

for iC=1:numCells
    bw=genCellRegion(cell_prop(iC), imageSize, cellRegion); %generate the mask for the current cell
    indSel=find(bw(ind));   %Note: bw(ind) will be the same length as fina_pts, with 0 or 1. So this returns the indices of points in the current cell
    pts_sel=final_pts(indSel,:);    %select out the points in the current cell
    cell_pts{iC}=pts_sel;
    posY=pts_sel(:,1)*pixSizeXY;    %Note: this converts the coord to nm
    posX=pts_sel(:,2)*pixSizeXY;
    posZ=pts_sel(:,3)*pixSizeZ;
    intensity=pts_sel(:,4);
    nSP=numel(posY);
    BGD=ones(nSP,1);
    residual=ones(nSP,1);
    if strcmp(p.method, 'GaussianFit')  %save the sigma's if the method is gaussian fit. 
        sigmaX=pts_sel(:,5)*pixSizeXY;  %this is in nm
        sigmaY=pts_sel(:,5)*pixSizeXY;
        sigmaZ=pts_sel(:,6)*pixSizeZ;
    else
        sigmaX=ones(nSP,1)*sigma_XY*pixSizeXY;
        sigmaY=ones(nSP,1)*sigma_XY*pixSizeXY;
        sigmaZ=ones(nSP,1)*sigma_Z*pixSizeZ;
    end
    %these are just to fill the blank required for FQ
    Cent_Y=ones(nSP,1);
    Cent_X=ones(nSP,1);
    Cent_Z=ones(nSP,1);
    MuY=ones(nSP,1);
    MuX=ones(nSP,1);
    MuZ=ones(nSP,1);
    ITERY_det=ones(nSP,1);
    Y_det=ceil(pts_sel(:,1));       %Note: all these are in pixel
    X_det=ceil(pts_sel(:,2));
    Z_det=ceil(pts_sel(:,3));
    Y_min=ceil(pts_sel(:,1))-boundSizeXY;
    Y_max=ceil(pts_sel(:,1))+boundSizeXY;
    Y_max(Y_min<1)=2*boundSizeXY+1;
    Y_min(Y_min<1)=1;    
    Y_min(Y_max>IMGSizeY)=IMGSizeY-2*boundSizeXY;
    Y_max(Y_max>IMGSizeY)=IMGSizeY;

    X_min=ceil(pts_sel(:,2))-boundSizeXY;
    X_max=ceil(pts_sel(:,2))+boundSizeXY;
    X_max(X_min<1)=2*boundSizeXY+1;
    X_min(X_min<1)=1;    
    X_min(X_max>IMGSizeX)=IMGSizeX-2*boundSizeXY;
    X_max(X_max>IMGSizeX)=IMGSizeX;
    
    Z_max=ceil(pts_sel(:,3))+boundSizeZ;
    Z_min=ceil(pts_sel(:,3))-boundSizeZ;
    Z_max(Z_min<1)=2*boundSizeZ+1;
    Z_min(Z_min<1)=1;
    Z_min(Z_max>IMGSizeZ)=IMGSizeZ-2*boundSizeZ;
    Z_max(Z_max>IMGSizeZ)=IMGSizeZ;
    INT_raw=ones(nSP,1);
    INT_filt=ones(nSP,1);
    SC_det=ones(nSP,1);
    SC_det_norm=ones(nSP,1);
    TH_det=ones(nSP,1);
    TH_fit=ones(nSP,1);
    spots_fit=[posY posX posZ intensity BGD residual sigmaX sigmaY sigmaZ ...
        Cent_Y Cent_X Cent_Z MuY MuX MuZ ITERY_det ];
    spots_det=[Y_det X_det Z_det ...
        Y_min Y_max X_min X_max Z_min Z_max INT_raw INT_filt SC_det ...
        SC_det_norm TH_det];
    thresh=struct('in',TH_fit,'all',TH_fit);
    cell_prop(iC).spots_fit=spots_fit;
    cell_prop(iC).spots_detected=spots_det;
    cell_prop(iC).thresh=thresh;
end
end
