function analysisTLS_RNA(filenames, p_RNA, maskType, cellRegion, verbose)

%the default parameter, could be changed
if nargin<3
    maskType=1;
end
if nargin<4
    cellRegion='cytosol';
end
if nargin<5
    verbose=false;
end

outlineFolder=filenames.outlineFolder;
resultFolder=filenames.resultFolder;
imageFolder=filenames.imageFolder;
RNAChannel=filenames.RNAChannel;
baseName=filenames.baseName;
fileNumber=filenames.fileNumber;

%% RNA detection
RNAFile=filenames.fileFunc(RNAChannel, baseName, fileNumber);
imageFile=fullfile(imageFolder,[RNAFile '.tif']);
outlineFile=fullfile(outlineFolder, [RNAFile,'__outline.txt']);
resultFile=fullfile(resultFolder,[RNAFile,'.txt']);

matResults=fullfile(resultFolder, [RNAFile, '.mat']);

%%read files
%Read in the RNA outline File
%[~,para_microscope,file_names,~,version]=FQ_load_results_WRAPPER_v2(outlineFile,'');
imf=imfinfo(imageFile);
nz=numel({imf.Height});
img=tiffread5(imageFile,1,nz);

%% Fitting the RNA
[cell_pts_RNA,~,cell_prop_RNA,final_pts_RNA,~, p_RNA]=uLocalizeFQ(img,p_RNA,outlineFile, ...
    resultFile, maskType, cellRegion, 'verbose', verbose);

if strcmpi(p_RNA.method,'GaussianFit')
    disp(['RNA channel: sigma_xy=', num2str(median(final_pts_RNA(:,5))), '; sigma_z=', num2str(median(final_pts_RNA(:,6)))]);
end
if verbose
    disp(['total RNAs: ', num2str(arrayfun(@(x) size(cell_prop_RNA(x).spots_fit,1), 1:numel(cell_prop_RNA)))]);
    disp(['single intensity: ', num2str(arrayfun(@(x) median(cell_prop_RNA(x).spots_fit(:,4)), 1:numel(cell_prop_RNA)),'%8.1f, ')]);
end

save(matResults,'cell_pts_RNA', 'cell_prop_RNA', 'final_pts_RNA', 'p_RNA');
end