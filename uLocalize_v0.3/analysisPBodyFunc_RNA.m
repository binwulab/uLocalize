function [cell_prop_RNA, cell_ts_RNA, cell_granule_RNA, p_RNA]=analysisPBodyFunc_RNA( ...
    filenames, p_RNA, showRNA, displayRange, calibration, maskType, cellRegion, verbose)
%the default parameter, could be changed
if nargin<3
    showRNA=false;
end
if nargin<4 || isempty(displayRange)
    displayRange=[0.1,0.2];
end
if nargin<5
    calibration=false;
end
if nargin<6
    maskType=1;
end
if nargin<7
    cellRegion='cytosol';
end
if nargin<8
    verbose=false;
end

outlineFolder=filenames.outlineFolder;
resultFolder=filenames.resultFolder;
imageFolder=filenames.imageFolder;
baseName=filenames.baseName;
fileNumber=filenames.fileNumber;
pBodyChannel=filenames.pBodyChannel;
mRNAChannel=filenames.mRNAChannel;

RNAFile=filenames.fileFunc(mRNAChannel, baseName, fileNumber);
imageFile=fullfile(imageFolder, [RNAFile, '.tif']);
outlineFile=fullfile(outlineFolder,[RNAFile,'_outline.txt']);
resultFile=fullfile(resultFolder,[RNAFile,'.txt']);

matResults=fullfile(resultFolder, [RNAFile, '.mat']);

%% Need to load pBody files for some information to generate RNA outline and quantification
pBodyFile=filenames.fileFunc(pBodyChannel, baseName, fileNumber);
pBody_matResults=fullfile(resultFolder,[pBodyFile,'.mat']);
outlineFile_pBody=fullfile(outlineFolder, [pBodyFile,'_outline.txt']);
[~,para_microscope,file_names,~,version]=FQ_load_results_WRAPPER_v2(outlineFile_pBody,''); %Need this info to generate outline
load(pBody_matResults, 'cell_granule_pBody', 'cell_prop_pBody', 'p_pBody');   %Need this info for detecting RNA granule

%% generate the mRNA outlines
% taking pBody position out and quantify separately. 
% To do that, input the pBody position as TS positions. For TS quantification, we should add a convergent criterion
cell_prop_RNA=cell_prop_pBody;
for i=1:numel(cell_prop_RNA)
    cc=cell_granule_pBody(i);   %the granule cc for i^th cell
    pos_TS=repmat(struct('label','','x',[],'y',[]), 1, cc.NumObjects);
    for j=1:cc.NumObjects
        %bndBox=getObjBndBox(cc.PixelIdxList{j}, cc.ImageSize,
        %p_RNA.bg_extension); %LBlake changed this to p_pBody.bg_extension
        bndBox=getObjBndBox(cc.PixelIdxList{j}, cc.ImageSize, p_pBody.bg_extension);
        xd=[bndBox(1,2), bndBox(2,2), bndBox(2,2), bndBox(1,2)];
        yd=[bndBox(1,1), bndBox(1,1), bndBox(2,1), bndBox(2,1)];
        pos_TS(j)=struct('label', num2str(j), 'x',xd, 'y', yd);
    end
    cell_prop_RNA(i).pos_TS=pos_TS;
    cell_prop_RNA(i).spots_fit=[];
    cell_prop_RNA(i).spots_det=[];
    cell_prop_RNA(i).thresh=[];
end
%save the RNA outline file
file_names.raw=[RNAFile, '.tif'];
RNAPara=struct('path_save', outlineFile, 'cell_prop', cell_prop_RNA, ...
    'par_microscope', para_microscope, 'file_names',file_names, 'version', version, 'flag_type', 'spot');
FQ_save_results_v1(outlineFile, RNAPara);
clear('pBodyFile', 'pBody_matResults', 'outlineFile_pBody', 'cell_prop_pBody', ...
    'pos_TS', 'xd','yd', 'cell_prop_RNA')

%% Quantification of RNA channel
%read in the RNA image file
imf=imfinfo(imageFile);
nz=numel({imf.Height});
img=mytiffread(imageFile,1:nz);

%%Fitting the RNA
[cell_pts_RNA,cell_ts_RNA,cell_prop_RNA,final_pts_RNA,~, p_RNA, cell_granule_RNA] = ...
    uLocalizeFQ(img,p_RNA,outlineFile,resultFile, maskType, cellRegion, ...
    'granuleCC', cell_granule_pBody, 'calibration', calibration,'verbose', verbose);
if showRNA 
    showRNAULocalize(img,  cell_prop_RNA, displayRange, cell_granule_RNA)
end
if strcmpi(p_RNA.method,'GaussianFit')
    disp(['RNA channel: sigma_xy=', num2str(median(final_pts_RNA(:,5))), '; sigma_z=', num2str(median(final_pts_RNA(:,6)))]);
end
if strcmpi(p_RNA.detectionMode,'CC') && calibration
    disp(['numDilation:',num2str(p_RNA.numDilation), ' bg_extension:', num2str(p_RNA.bg_extension), ...
        ' ISum2IFitRatio: ', num2str(p_RNA.ISum2IFitRatio), ' aMax:', num2str(p_RNA.aMax)]);
end
if verbose
    if strcmpi(p_RNA.detectionMode,'CC') 
        disp(['total RNAs: ', num2str(arrayfun(@(x) sum(cell_prop_RNA(x).single.fitStatus) + sum(cell_prop_RNA(x).multiple.nRNA) + sum(cell_prop_RNA(x).spotBG.nRNA), 1:numel(cell_prop_RNA)))]);    
        disp(['single intensity: ', num2str(arrayfun(@(x) median(cell_prop_RNA(x).single.fitRes(:,4)), 1:numel(cell_prop_RNA)),'%8.1f, ')]);
    else
        disp(['total RNAs: ', num2str(arrayfun(@(x) size(cell_prop_RNA(x).spots_fit,1), 1:numel(cell_prop_RNA)))]);
        disp(['single intensity: ', num2str(arrayfun(@(x) median(cell_prop_RNA(x).spots_fit(:,4)), 1:numel(cell_prop_RNA)),'%8.1f, ')]);
    end
end
save(matResults,'cell_pts_RNA', 'cell_ts_RNA', 'cell_prop_RNA', 'final_pts_RNA', 'cell_granule_RNA', 'p_RNA');
end