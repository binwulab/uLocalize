function analysisPBodyFunc_pBody(filenames, p_pBody, maskType, cellRegion, verbose)

if nargin<3
    maskType=1;
end
if nargin<4
    cellRegion='cytosol';
end
if nargin<5
    verbose=false;
end
%the default parameter, could be changed

outlineFolder=filenames.outlineFolder;
resultFolder=filenames.resultFolder;
imageFolder=filenames.imageFolder;
pBodyChannel=filenames.pBodyChannel;
baseName=filenames.baseName;
fileNumber=filenames.fileNumber;

%% p-body detection
pBodyFile=filenames.fileFunc(pBodyChannel, baseName, fileNumber);
imageFile=fullfile(imageFolder,[pBodyFile '.tif']);
outlineFile=fullfile(outlineFolder, [pBodyFile,'_outline.txt']);
resultFile=fullfile(resultFolder,[pBodyFile,'.txt']);
%filteredImage=fullfile(imageFolder,[pBodyFile,'_filtered.mat']);   %Not to save the filtered image

%matResults=fullfile(resultFolder, [baseName, num2str(fileNumber), '.mat']);
matResults=fullfile(resultFolder, [pBodyFile, '.mat']);

%%read files
%Read in the pBody File
% [cell_prop_pBody,para_microscope,file_names,~,version]=FQ_load_results_WRAPPER_v2(outlineFile,'');
imf=imfinfo(imageFile);
nz=numel({imf.Height});
img=tiffread5(imageFile,1,nz);

%%fitting the pBody
[cell_pts_pBody, ~, cell_prop_pBody, final_pts_pBody, ~, p_pBody, cell_granule_pBody]= ...
    uLocalizeFQ(img,p_pBody,outlineFile,resultFile, maskType,cellRegion, 'verbose', verbose);

save(matResults,'cell_pts_pBody', 'cell_prop_pBody', 'final_pts_pBody', 'cell_granule_pBody', 'p_pBody');
end


