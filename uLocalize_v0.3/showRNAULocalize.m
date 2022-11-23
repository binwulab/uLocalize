function showRNAULocalize(img,  cell_prop_RNA, displayRange, cell_granule_RNA)

if nargin<3
    displayRange=[0.05,0.2];
end
imgMP=max(img,[],3);
mx=double(max(imgMP(:)));
mn=double(min(imgMP(:)));
imshow(imgMP,[mn+(mx-mn)*displayRange(1), mn+(mx-mn)*displayRange(2)]);
hold on;
for i=1:numel(cell_prop_RNA)
    if isfield(cell_prop_RNA,'multiple')    % CC mode: strcmpi(p_RNA.detectionMode,'CC')
        plot(cell_prop_RNA(i).single.fitRes(:,2), cell_prop_RNA(i).single.fitRes(:,1), 'ro');
        plot(cell_prop_RNA(i).multiple.fitRes(:,2), cell_prop_RNA(i).multiple.fitRes(:,1), 'yo');
        plot(cell_prop_RNA(i).spotBG.fitRes(:,2), cell_prop_RNA(i).spotBG.fitRes(:,1), 'co');
    else
        plot(cell_prop_RNA(i).spots_detected(:,2), cell_prop_RNA(i).spots_detected(:,1), 'ro');
    end
    if ~isempty(cell_granule_RNA)
        plot(cell_granule_RNA(i).ObjCtr(:,2), cell_granule_RNA(i).ObjCtr(:,1), 'gs', 'MarkerSize', 10);
    end
end
hold off
end