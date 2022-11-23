function analysisPBodyFuncV2(filenames, p_pBody, p_RNA, showRNA, displayRange, calibration)

%the default parameter, could be changed
maskType=1;
cellRegion='cytosol';
forceFilter=true;
if nargin<5
    showRNA=false;
end
if nargin<6 || isempty(displayRange)
    displayRange=[0.1,0.2];
end
if nargin<7
    calibration=false;
end

rootFolder=filenames.rootFolder;
outlineFolder=filenames.outlineFolder;
resultFolder=filenames.resultFolder;
imageFolder=filenames.imageFolder;
pBodyChannel=filenames.pBodyChannel;
mRNAChannel=filenames.mRNAChannel;
baseName=filenames.baseName;
fileNumber=filenames.fileNumber;

%% p-body detection
pBodyFile=filenames.fileFunc(pBodyChannel, baseName, fileNumber);
%pBodyFile=[pBodyChannel, baseName, num2str(fileNumber, '%02d')];
imageFile=fullfile(imageFolder,[pBodyFile '.tif']);
outlineFile=fullfile(outlineFolder, [pBodyFile,'_outline.txt']);
resultFile=fullfile(resultFolder,[pBodyFile,'.txt']);
filteredImage=fullfile(imageFolder,[pBodyFile,'_filtered.mat']);

matResults=fullfile(resultFolder, [baseName, num2str(fileNumber), '.mat']);

%%read files
%Read in the pBOdy File
[cell_prop_pBody,para_microscope,file_names,~,version]=FQ_load_results_WRAPPER_v2(outlineFile,'');
imf=imfinfo(imageFile);
nz=numel({imf.Height});
img=tiffread5(imageFile,1,nz);
%read in the filtered image
% if exist(filteredImage,'file')==2 && ~forceFilter
%     load(filteredImage);    %It contained the variable smooth, which is the smoothed image
% else
%     smooth=[];
% end

%%fitting the pBody
% thresh=thresh_pBody;
% aMin=aMin_pBody;
% eMax=eMax_pBody;
% bg_ext=bg_extension_pBody;
% filterSigma=filterSigma_pBody;
% 
% p=uLocalizeInitPara('numdim',3, 'aMin', aMin, 'eMax', eMax, 'detectionMode', 'granule', 'thresh', thresh, ...
%     'filterMethod', 'LOGRAJ', 'filterSigma', filterSigma, 'thickness',1, 'bg_extension', bg_ext);
[cell_pts_pBody, ~, cell_prop_pBody, final_pts_pBody, ~, p_pBody, cell_granule_pBody]=uLocalizeFQ(img,p_pBody,outlineFile,resultFile, maskType, ...
    cellRegion);
%[final_pts, p, ~, cc]=uLocalizeDetection(img, p, maskType, mask, 'thresh', thresh);
% if exist(filteredImage,'file')~=2 || forceFilter
%     save(filteredImage,'smooth');
% end

%% mRNA detection
% RNAFile=[mRNAChannel, baseName, num2str(fileNumber, '%02d')];
RNAFile=filenames.fileFunc(mRNAChannel, baseName, fileNumber);
imageFile=fullfile(imageFolder, [RNAFile, '.tif']);
outlineFile=fullfile(outlineFolder,[RNAFile,'_outline.txt']);
resultFile=fullfile(resultFolder,[RNAFile,'.txt']);
filteredImage=fullfile(imageFolder,[RNAFile,'_filtered.mat']);

%read in the RNA image file
imf=imfinfo(imageFile);
nz=numel({imf.Height});
%img=mytiffread(imageFile,1,nz);
img=mytiffread(imageFile,1:nz);
forceFilter=false;
%read in the filtered image
% if exist(filteredImage,'file')==2 && ~forceFilter
%     load(filteredImage);    %It contained the variable smooth, which is the smoothed image
% else
%     smooth=[];
% end

%%Making the RNA channel outline, taking pBody position out and
%%quantify separately. To do that, input the pBody position as TS positions. For TS quantification, we should add a convergent criterion
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

%%Fitting the RNA
maskTypeRNA=1; % for RNA only. use set value for std LBlake
[cell_pts_RNA,cell_ts_RNA,cell_prop_RNA,final_pts_RNA,~, p_RNA, cell_granule_RNA]=uLocalizeFQ(img,p_RNA,outlineFile,resultFile, maskTypeRNA, ...
    cellRegion, 'granuleCC', cell_granule_pBody, 'calibration', calibration);
% save filtered image, so no need to calculate it again
% if exist(filteredImage,'file')~=2 || forceFilter
%     save(filteredImage,'smooth');
% end
if showRNA 
    imgMP=max(img,[],3);
    mx=double(max(imgMP(:)));
    mn=double(min(imgMP(:)));
    imshow(imgMP,[mn+(mx-mn)*displayRange(1), mn+(mx-mn)*displayRange(2)]);
    hold on;
    for i=1:numel(cell_prop_RNA)
        if strcmpi(p_RNA.detectionMode,'CC')
            plot(cell_prop_RNA(i).single.fitRes(:,2), cell_prop_RNA(i).single.fitRes(:,1), 'ro');
            plot(cell_prop_RNA(i).multiple.fitRes(:,2), cell_prop_RNA(i).multiple.fitRes(:,1), 'yo');
        else
            plot(cell_pts_RNA{i}(:,2), cell_pts_RNA{i}(:,1), 'ro');
        end
        plot(cell_granule_RNA(i).ObjCtr(:,2), cell_granule_RNA(i).ObjCtr(:,1), 'gs', 'MarkerSize', 10);
    end
    hold off
end
if strcmpi(p_RNA.method,'GaussianFit')
    disp(['RNA channel: sigma_xy=', num2str(median(final_pts_RNA(:,5))), '; sigma_z=', num2str(median(final_pts_RNA(:,6)))]);
end
if strcmpi(p_RNA.detectionMode,'CC') && calibration
    disp(['numDilation:',num2str(p_RNA.numDilation), ' bg_extension:', num2str(p_RNA.bg_extension), ' ISum2IFitRatio: ', num2str(p_RNA.ISum2IFitRatio)]);
end
save(matResults,'cell_pts_pBody', 'cell_prop_pBody', 'final_pts_pBody', 'cell_granule_pBody', 'p_pBody', 'cell_pts_RNA', 'cell_ts_RNA', 'cell_prop_RNA', 'final_pts_RNA', 'cell_granule_RNA', 'p_RNA');

end