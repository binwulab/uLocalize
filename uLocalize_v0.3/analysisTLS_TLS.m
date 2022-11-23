function [cell_prop_Protein, cell_ts_Protein, p_Protein]=analysisTLS_TLS( ...
    filenames, p_Protein,  maskType, cellRegion, verbose)
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
baseName=filenames.baseName;
fileNumber=filenames.fileNumber;
ProteinChannel=filenames.ProteinChannel;
RNAChannel=filenames.RNAChannel;


%% Need to load RNA files for some information to generate Protein outline and quantification
RNAFile=filenames.fileFunc(RNAChannel, baseName, fileNumber);
RNA_matResults=fullfile(resultFolder,[RNAFile,'.mat']);
outlineFile_RNA=fullfile(outlineFolder, [RNAFile,'__outline.txt']);
[~,para_microscope,file_names,~,version]=FQ_load_results_WRAPPER_v2(outlineFile_RNA,''); %Need this info to generate outline
load(RNA_matResults, 'cell_pts_RNA', 'cell_prop_RNA');   %Need this info for detecting RNA granule

%% Protein detection
ProteinFile=filenames.fileFunc(ProteinChannel, baseName, fileNumber);
imageFile=fullfile(imageFolder, [ProteinFile, '.tif']);
outlineFile=fullfile(outlineFolder,[ProteinFile,'__outline.txt']);
resultFile=fullfile(resultFolder,[ProteinFile,'.txt']);

matResults=fullfile(resultFolder, [ProteinFile, '.mat']);

%% generate the Protein outlines
% taking RNA position out and quantify separately. 
% To do that, input the RNA position as TS positions. For TS quantification, we should add a convergent criterion
cell_prop_Protein=cell_prop_RNA;
TS_cutwidth=p_Protein.TS_cutwidth(1);
for i=1:numel(cell_prop_Protein)
    pts=cell_pts_RNA{i};
    if ~isempty(pts)
        yd=round(pts(:,1));
        xd=round(pts(:,2));
        nd=numel(xd);
        %prepare the box surrounding the TS and save it in a cell array
        xd=[xd-TS_cutwidth, xd+TS_cutwidth, xd+TS_cutwidth, xd-TS_cutwidth];
        yd=[yd-TS_cutwidth, yd-TS_cutwidth, yd+TS_cutwidth, yd+TS_cutwidth];
        xd=mat2cell(xd, ones(1,nd));
        yd=mat2cell(yd, ones(1,nd));
        pos_TS=struct('label',mat2cell(num2str((1:nd)'),ones(1,nd)),'x',xd,'y',yd); %distribute the box into TS.
        cell_prop_Protein(i).pos_TS=pos_TS';
    end
    cell_prop_Protein(i).spots_fit=[];  %clear the structure
    cell_prop_Protein(i).spots_det=[];
    cell_prop_Protein(i).thresh=[];
end
%save the Protein outline file
file_names.raw=[ProteinFile, '.tif'];
ProteinPara=struct('path_save', outlineFile, 'cell_prop', cell_prop_Protein, ...
    'par_microscope', para_microscope, 'file_names',file_names, 'version', version, 'flag_type', 'spot');
FQ_save_results_v1(outlineFile, ProteinPara);
clear('RNAFile', 'RNA_matResults', 'outlineFile_RNA', 'cell_prop_RNA', ...
    'pos_TS', 'xd','yd')

%% Quantification of RNA channel
%read in the RNA image file
imf=imfinfo(imageFile);
nz=numel({imf.Height});
img=mytiffread(imageFile,1:nz);

%%Fitting the RNA
[cell_pts_Protein,cell_ts_Protein,cell_prop_Protein,final_pts_Protein,~, p_Protein] = ...
    uLocalizeFQ(img,p_Protein,outlineFile,resultFile, maskType, cellRegion, 'verbose', verbose);

if strcmpi(p_Protein.method,'GaussianFit')
    disp(['RNA channel: sigma_xy=', num2str(median(final_pts_Protein(:,5))), '; sigma_z=', num2str(median(final_pts_Protein(:,6)))]);
end
if verbose
    disp(['total Proteins: ', num2str(arrayfun(@(x) size(cell_prop_Protein(x).spots_fit,1), 1:numel(cell_prop_Protein)))]);
    disp(['single intensity: ', num2str(arrayfun(@(x) median(cell_prop_Protein(x).spots_fit(:,4)), 1:numel(cell_prop_Protein)),'%8.1f, ')]);
end
save(matResults,'cell_pts_Protein', 'cell_ts_Protein', 'cell_prop_Protein', 'final_pts_Protein', 'p_Protein');
end