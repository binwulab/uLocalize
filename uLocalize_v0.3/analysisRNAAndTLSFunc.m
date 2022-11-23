function analysisRNAAndTLSFunc(filenames, fileNumber, p_RNA, p_Protein, cutsize_TS)

%the default parameter, could be changed
maskType=1;
cellRegion='cytosol';
forceFilter=false;
filterMethod='LOGRAJ';
if nargin < 5
    cutsize_TS=5;
end

rootFolder=filenames.rootFolder;
outlineFolder=filenames.outlineFolder;
resultFolder=filenames.resultFolder;
imageFolder=filenames.imageFolder;
ProteinChannel=filenames.ProteinChannel;
RNAChannel=filenames.RNAChannel;
baseName=filenames.baseName;

%% RNA detection
RNAFile=filenames.fileFunc(RNAChannel, baseName, fileNumber);
%pBodyFile=[pBodyChannel, baseName, num2str(fileNumber, '%02d')];
imageFile=fullfile(imageFolder,[RNAFile '.tif']);
outlineFile=fullfile(outlineFolder, [RNAFile,'__outline.txt']);
resultFile=fullfile(resultFolder,[RNAFile,'.txt']);
filteredImage=fullfile(imageFolder,[RNAFile,'_filtered.mat']);

matResults=fullfile(resultFolder, [baseName, num2str(fileNumber), '.mat']);

%%read files
%Read in the pBOdy File
[~,para_microscope,file_names,~,version]=FQ_load_results_WRAPPER_v2(outlineFile,'');
imf=imfinfo(imageFile);
nz=numel({imf.Height});
img=tiffread5(imageFile,1,nz);
%read in the filtered image
if exist(filteredImage,'file')==2 && ~forceFilter
    load(filteredImage);    %It contained the variable smooth, which is the smoothed image
else
    smooth=[];
end
[cell_pts_RNA,~,cell_prop_RNA,final_pts_RNA,smooth]=uLocalizeFQ(img,p_RNA,outlineFile,resultFile, maskType, ...
    cellRegion,'filteredImage',smooth,'forceFilter',forceFilter);
if exist(filteredImage,'file')~=2 || forceFilter
    save(filteredImage,'smooth');
end
if strcmpi(p_RNA.method,'GaussianFit')
    disp(['RNA channel: sigma_xy=', num2str(median(final_pts_RNA(:,5))), '; sigma_z=', num2str(median(final_pts_RNA(:,6)))]);
end


    %% defining protein channel
    ProteinFile=filenames.fileFunc(ProteinChannel, baseName, fileNumber);
    imageFile=fullfile(imageFolder, [ProteinFile, '.tif']);
    outlineFile=fullfile(outlineFolder,[ProteinFile,'__outline.txt']);
    resultFile=fullfile(resultFolder,[ProteinFile,'.txt']);
    filteredImage=fullfile(imageFolder,[ProteinFile,'_filtered.mat']);

    %% Making protein outline, taking RNA position as candidate translation sites
    cell_prop_Protein=cell_prop_RNA;
    for i=1:numel(cell_prop_Protein)
        pts=cell_pts_RNA{i};
        if ~isempty(pts)
            yd=round(pts(:,1));
            xd=round(pts(:,2));
            nd=numel(xd);
            %prepare the box surrounding the TS and save it in a cell array
            xd=[xd-cutsize_TS, xd+cutsize_TS, xd+cutsize_TS, xd-cutsize_TS];
            yd=[yd-cutsize_TS, yd-cutsize_TS, yd+cutsize_TS, yd+cutsize_TS];
            xd=mat2cell(xd, ones(1,nd));
            yd=mat2cell(yd, ones(1,nd));
            pos_TS=struct('label',mat2cell(num2str((1:nd)'),ones(1,nd)),'x',xd,'y',yd); %distribute the box into TS. 
            cell_prop_Protein(i).pos_TS=pos_TS';
        end
        cell_prop_Protein(i).spots_fit=[];  %clear the structure
        cell_prop_Protein(i).spots_det=[];
        cell_prop_Protein(i).thresh=[];
    end
    % save the outline file
    file_names.raw=[ProteinFile, '.tif'];
    proteinPara=struct('path_save', outlineFile, 'cell_prop', cell_prop_Protein, ...
        'par_microscope', para_microscope, 'file_names',file_names, 'version', version, 'flag_type', 'spot');
    FQ_save_results_v1(outlineFile, proteinPara);

    %% Read in the protein image
    imf=imfinfo(imageFile);
    nz=numel({imf.Height});
    img=tiffread5(imageFile,1,nz);
    forceFilter=false;
    %read in the filtered image
    if exist(filteredImage,'file')==2 && ~forceFilter
        load(filteredImage);    %It contained the variable smooth, which is the smoothed image
    else
        smooth=[];
    end

    %% fitting the protein image
    [cell_pts_Protein,cell_ts_Protein,cell_prop_Protein,final_pts_Protein,smooth]=uLocalizeFQ(img,p_Protein,outlineFile,resultFile, maskType, ...
        cellRegion,'filteredImage',smooth,'forceFilter',forceFilter);
    % save filtered image, so no need to calculate it again
    if exist(filteredImage,'file')~=2 || forceFilter
        save(filteredImage,'smooth');
    end
    if strcmpi(p_Protein.method,'GaussianFit')
        disp(['Protein channel: sigma_xy=', num2str(median(final_pts_Protein(:,5))), '; sigma_z=', num2str(median(final_pts_Protein(:,6)))]);
    end
    save(matResults,'cell_pts_RNA', 'cell_prop_RNA', 'final_pts_RNA', 'cell_pts_Protein', 'cell_ts_Protein', 'cell_prop_Protein', 'final_pts_Protein'); 
end