% rootFolder='C:\Users\Bin\OneDrive - Johns Hopkins University\data\Test\FISH\';
rootFolder='E:\Data\temp\Gosia\test\';
imageFolder=[rootFolder];
outlineFolder=[rootFolder 'Outlines'];
resultFolder=[rootFolder 'Results'];
imgFile=fullfile(imageFolder,'C1-A60_NT_1.tif');
outlineFile=fullfile(outlineFolder,'C1-A60_NT_1__outline.txt');
filteredImage=fullfile(imageFolder,'C1-A60_NT_1_filtered.mat');
outputFile=fullfile(resultFolder,'C1-A60_NT_1.txt');
imf=imfinfo(imgFile);
nz=numel({imf.Height});
img=tiffread5(imgFile,1,nz);
%%
p=uLocalizeInitPara('numdim',3, 'sigma_xy', 1.5, 'sigma_z',1.8,'thresh',20, 'filterMethod', 'LOGFFT', 'filterSigma',1.4, 'cutsize',3, 'thickness',2,'method','GaussianMask');
maskType=1;
forceFilter=false;
if exist(filteredImage,'file')==2 && ~forceFilter
    load(filteredImage);    %It contained the variable smooth, which is the smoothed image
else
    smooth=[];
end
p.denseMode=1;
[cell_pts,cell_ts,cell_prop,final_pts,smooth]=uLocalizeFQ(img,p,outlineFile,outputFile, maskType,'cytosol','filteredImage',smooth,'forceFilter',forceFilter,'calibration',true);
if exist(filteredImage,'file')~=2 || forceFilter
    save(filteredImage,'smooth');
end

%% to calibrate the sigma
maskType=1;
[cell_pts,cell_ts,cell_prop,final_pts,smooth]=uLocalizeFQ(img,p,outlineFile,outputFile, maskType,'cytosol','filteredImage', smooth,'thresh',13, 'method','GaussianFit');

%%
imshow(max(img,[],3),[200,2000]);
hold on;
plot(final_pts(:,2),final_pts(:,1),'sq');