resultDir='Z:\users\nliving5\2020\E4.80 FISH-IF UTR and ORF Cell Lines\';

for image=[2:16]
mRNAResult=fullfile(resultDir,'Results\ST\', ['C1-ST_' num2str(image) '__outline_spots_201217.txt']);
[cell_prop_RNA, par_microscope_RNA, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(mRNAResult,'');
pixelRNA=par_microscope_RNA.pixel_size.xy;
proteinOutline=fullfile(resultDir,'Outlines\', ['C3-ST_' num2str(image) '__outline.txt']);

[cell_prop_Protein, par_microscope_Protein, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(proteinOutline,'');
crop_xy=5;


for i=1:numel(cell_prop_Protein)
    if ~isempty(cell_prop_RNA(i).spots_fit)
    yd=round(cell_prop_RNA(i).spots_fit(:,1)/pixelRNA-0.5)+1;
    xd=round(cell_prop_RNA(i).spots_fit(:,2)/pixelRNA-0.5)+1;
    pos_TS=repmat(struct('label','','x',[],'y',[],'z',[]),[numel(xd),1]);
    for j=1:numel(xd)
        pos_TS(j).label=['TxS_', num2str(j)];
        pos_TS(j).x=[xd(j)-crop_xy, xd(j)+crop_xy, xd(j)+crop_xy, xd(j)-crop_xy];
        pos_TS(j).y=[yd(j)-crop_xy, yd(j)-crop_xy, yd(j)+crop_xy, yd(j)+crop_xy];
    end
    cell_prop_Protein(i).pos_TS=pos_TS;
    end
end

proteinOutlineTS=fullfile(resultDir,'Outlines\', ['C3-ST_' num2str(image) '_TS_outline.txt']);
proteinPara=struct('path_save', proteinOutlineTS, 'cell_prop', cell_prop_Protein, 'par_microscope', par_microscope_Protein, 'file_names', file_names, 'version', version, 'flag_type', '');

FQ_save_results_v1(proteinOutlineTS, proteinPara);
end