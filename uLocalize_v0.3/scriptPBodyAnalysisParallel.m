%scriptPBodyAnalysisParallel

%% p-Body detection parameters
maskType=1;
cellRegion='cytosol';
sigma_xy=1.44;
sigma_z=1.4;
cutsize=3;
thickness=1;
% threshUnit='sd';
threshUnit='absolute';
threshLevel=70;
aMin=7;
eMax=0.9;
bg_extension=[3,2];
filterSigma_pBody=1.*sigma_xy;
filterWindowSize = max(2*round(2*filterSigma_pBody)+1,5);
bgCorrMethod = 'quad';
p_pBody=uLocalizeInitPara('numdim',3, 'cutsize', cutsize, 'thickness',thickness, 'threshUnit', threshUnit, 'threshLevel', threshLevel, ...
    'sigma_xy', sigma_xy, 'sigma_z', sigma_z, 'filterMethod', 'LOGRAJ', 'filterSigma', filterSigma_pBody,  'filterWindowSize', filterWindowSize, ...
    'aMin', aMin, 'eMax', eMax, 'detectionMode', 'granule', 'bg_extension', bg_extension, 'bgCorrMethod', bgCorrMethod);
%% RNA detection parameters
% sigma_xy=1.21;
% sigma_z=0.95;

%cy3
sigma_xy=1.03;
sigma_z=0.96;

%cy5
% sigma_xy=1.19;
% sigma_z=1.16;

% for 2 hr rapa drb: sigma_xy=1.2801; sigma_z=1.5171

%RNA channel: sigma_xy=1.0608; sigma_z=0.75711
%RNA channel: sigma_xy=1.066; sigma_z=0.8297
%RNA channel: sigma_xy=1.0795; sigma_z=1.1288

% threshUnit='sd';
threshUnit='absolute';
threshLevel=90;
cutsize=2;
thickness=1;
detectionMode='CC';
aMin=5;
aMax=30;
eMax=0.8;
numDilation=3;
bg_extension=[2, 1];
ISum2IFitRatio=0.9;
fitmethod='GaussianMask';
% fitmethod='GaussianFit';   %use this method to calibrate the sigma
filterMethod='LOGRAJ';
filterSigma_RNA=sigma_xy;
filterWindowSize= max(2*round(2*filterSigma_RNA)+1,5);
NormByCellLimit=10;
bgCorrMethod = 'plane';
calibration=false;
showRNA=false;
displayRange=[0.05,0.2];
verbose = true;
InternalParallel = false;

p_RNA=uLocalizeInitPara('numdim',3, 'cutsize',cutsize, 'thickness',thickness, 'threshUnit', threshUnit, 'threshLevel', threshLevel, ...
    'sigma_xy', sigma_xy, 'sigma_z',sigma_z,  'filterMethod',filterMethod,'filterSigma', filterSigma_RNA, 'filterWindowSize', filterWindowSize,...
    'InternalParallel',InternalParallel, 'detectionMode', detectionMode, 'method',fitmethod, 'NormByCellLimit', NormByCellLimit, 'bgCorrMethod', bgCorrMethod, ...
    'aMin', aMin, 'aMax', aMax, 'eMax', eMax, ...
    'numDilation', numDilation, 'bg_extension', bg_extension, 'ISum2IFitRatio', ISum2IFitRatio);

rootFolder='\\10.17.0.34\BinWuLab2\users\lblake8\BinToRun\210622\';
outlineFolder=[rootFolder 'Outlines\'];
resultFolder=[rootFolder '\Results\test\'];
imageFolder=rootFolder;
pBodyChannel='_AF750';
mRNAChannel='_Cy3';

filenamesAll=defineFileNames(rootFolder, outlineFolder, resultFolder, imageFolder, pBodyChannel, mRNAChannel);
% parpool('local',32);

%% for debugging purpose
% analysisPBodyFunc_pBody(filenamesAll(1), p_pBody, maskType, cellRegion, verbose);
% analysisPBodyFunc_RNA(filenamesAll(1), p_RNA, showRNA, displayRange, calibration, maskType, cellRegion, verbose);
% return;

errFiles=[];
errFileNames={};

%% p_body detection
for i=1:numel(filenamesAll)  %This doesn't have to be executed every time, can be commented out once p_body detection is satisfactory
    try
        analysisPBodyFunc_pBody(filenamesAll(i), p_pBody, maskType, cellRegion, verbose);
    catch ME
        disp(['ERROR: ' filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]);
        errFiles=[errFiles; filenamesAll(i)];
        errFileNames=[errFileNames, [filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]];
        continue
    end
end

%% RNA detection

% for i=1:numel(filenamesAll)
%     try
%         [cp,ct,cg,p]=analysisPBodyFunc_RNA(filenamesAll(i), p_RNA, showRNA, displayRange, calibration, maskType, cellRegion, verbose);
%     catch ME
%         disp(['ERROR: ' filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]);
%         errFiles=[errFiles; filenamesAll(i)];
%         errFileNames=[errFileNames, [filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]];
%         continue
%     end
% end

save(fullfile(rootFolder, 'errFiles.mat'), 'errFiles','errFileNames');
% delete(gcp);

function filenamesAll=defineFileNames(rootFolder, outlineFolder, resultFolder, imageFolder, pBodyChannel, mRNAChannel)
filenamesAll=[];
%add all files here

% baseName='mbsmef_1_platea_12hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_2_platea_9hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_3_platea_6hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_4_platea_3hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% % fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% fileNumbers={1 2 3};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_5_platea_2hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={2 3 4};
% %fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% % % 
% baseName='mbsmef_6_platea_1hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={2};
% % fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_7_platea_30minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% %fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% fileNumbers={2};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_8_platea_15minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% % fileNumbers={1 2 3 4 5 6 7 8 9 10};
% fileNumbers={4};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_9_platea_10minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
baseName='mbsmef_10_platea_5minrapadrb_xy';
fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
fileNumbers={3};
filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
    'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
    'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_11_platea_steadystate_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% %fileNumbers={1 2 3 4 5 6 7 8 9 10};
% fileNumbers={2};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];

% baseName='mbsmef_12_platea_12hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% %fileNumbers={1};
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_13_platea_9hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_14_platea_6hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_15_platea_3hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_16_platea_2hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_17_platea_1hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_18_platea_30mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_19_platea_15mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_20_platea_10mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_21_platea_5mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_22_platea_12hractbrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_23_platea_9hractbrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_24_platea_6hractbrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_25_platea_2hractbrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_26_platea_12hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12 13};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_27_platea_9hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_28_platea_6hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_29_platea_2hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_30_platea_1hrdmso_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_31_platea_30mindmso_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_32_platea_15mindmso_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_33_platea_5mindmso_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% %fileNumbers={1 2 3 4 5 6 7 8 9 10 11 12};
% fileNumbers={1};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];

% baseName='mbsmef_34_plateb_12hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_35_plateb_9hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_36_plateb_6hrrapadrb002_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_37_plateb_3hrrapadrb001_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_38_plateb_2hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_39_plateb_1hrrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_40_plateb_30minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_41_plateb_15minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_42_plateb_10minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={2 3 4, 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_43_plateb_5minrapadrb_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_44_plateb_steadystate_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_45_plateb_12hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_46_plateb_9hrdrbonly001_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_47_plateb_6hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_48_plateb_3hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_49_plateb_2hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_50_plateb_1hrdrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_51_plateb_30mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_53_plateb_10mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_54_plateb_5mindrbonly_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_55_plateb_12hractbrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2, 4 5 6, 8 9, 11 12 13};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_56_plateb_9hractbrna_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_57_plateb_6hractbrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_58_plateb_2hractbrna_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3, 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_59_plateb_12hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5, 8, 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_60_plateb_9hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_61_plateb_6hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_62_plateb_2hrncrnai_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_63_plateb_1hrvehicle_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_64_plateb_30minvehiclei_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_65_plateb_15minvehicle_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9 10 11};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
% 
% baseName='mbsmef_66_plateb_5minvehicle_xy';
% fileFunc=@(channel, base, num) [base,num2str(num,'%02d'),channel]; %this needs to be modified according to the naming convention
% fileNumbers={1 2 3 4 5 6 7 8 9};
% filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
%     'resultFolder', resultFolder, 'imageFolder', imageFolder, 'pBodyChannel', pBodyChannel, ...
%     'mRNAChannel', mRNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
% filenamesAll=[filenamesAll, filenames];
end