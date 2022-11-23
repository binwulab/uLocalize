%scriptTLSAnalysisParallel

%% RNA detection parameters
maskType=1; %
cellRegion='cytosol';
sigma_xy=1.1718;
sigma_z=1.0303;
cutsize=2;
thickness=1;
threshUnit='sd'; threshLevel=11;
%threshUnit='absolute'; threshLevel=150;

detectionMode='LM';
fitmethod='GaussianMask';
%fitmethod='GaussianFit';   %use this method if want to calibrate the sigma

filterMethod='LOGRAJ';
filterSigma=sigma_xy;
filterWindowSize = max(2*round(2*filterSigma)+1,5);

bgCorrMethod = 'plane';
p_RNA=uLocalizeInitPara('numdim',3, 'cutsize', cutsize, 'thickness',thickness, 'threshUnit', threshUnit, 'threshLevel', threshLevel, ...
    'sigma_xy', sigma_xy, 'sigma_z', sigma_z, 'filterMethod', filterMethod, 'filterSigma', filterSigma,  'filterWindowSize', filterWindowSize, ...
    'detectionMode', detectionMode, 'method', fitmethod, 'bgCorrMethod', bgCorrMethod);
%% Protein detection parameters
sigma_xy=0.93578;
sigma_z=0.99501;
%thresh=8;
cutsize=3;
thickness=1;
TS_cutsize=5;   %The real cut_width is going to be ceil(TS_cutsize*sigma_xy)
% threshUnit='sd'; threshLevel=8;
threshUnit='absolute'; threshLevel=125;

detectionMode='LM';
fitmethod='GaussianMask';
%fitmethod='GaussianFit';   %use this method to calibrate the sigma

filterMethod='LOGRAJ';
filterSigma=sigma_xy;
filterWindowSize= max(2*round(2*filterSigma)+1,5);
bgCorrMethod = 'plane';
verbose = true;

p_Protein=uLocalizeInitPara('numdim',3, 'cutsize',cutsize, 'TS_cutsize', TS_cutsize, 'thickness',thickness, 'threshUnit', threshUnit, 'threshLevel', threshLevel, ...
    'sigma_xy', sigma_xy, 'sigma_z',sigma_z,  'filterMethod',filterMethod,'filterSigma', filterSigma, 'filterWindowSize', filterWindowSize,...
    'detectionMode', detectionMode, 'method',fitmethod, 'bgCorrMethod', bgCorrMethod);

%% define the folders and file names
rootFolder='Z:\users\nliving5\2021\E5.1 FISH-IF Repeat UTRs and ORFs\';
outlineFolder=[rootFolder 'Outlines'];
resultFolder=[rootFolder 'Results'];
imageFolder=[rootFolder, 'Images'];
ProteinChannel='C3-';
RNAChannel='C1-';

filenamesAll=defineFileNames(rootFolder, outlineFolder, resultFolder, imageFolder, ProteinChannel, RNAChannel);
% parpool('local',32);

%% for debugging purpose
% analysisTLS_RNA(filenamesAll(1), p_RNA, maskType, cellRegion, verbose);
% [cp,ct,p]=analysisTLS_TLS(filenamesAll(1), p_Protein, maskType, cellRegion, verbose);
% return;

errFiles=[];
errFileNames={};
%%  RNA detection, can change to parfor for parallel

parpool('local', 2); %make sure to run this to open specific number of threads, it can be run separately
parfor i=1:numel(filenamesAll)  %This doesn't have to be executed every time, can be commented out once RNA detection is satisfactory
    try
        analysisTLS_RNA(filenamesAll(i), p_RNA, maskType, cellRegion, verbose);
    catch ME
        disp(['ERROR: ' filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]);
        errFiles=[errFiles; filenamesAll(i)];
        errFileNames=[errFileNames, [filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]];
        continue
    end
end

delete(gcp);
%% Protein detection, can change to parfor for parallel

parpool('local', 2); %make sure to run this to open specific number of threads, it can be run separately
parfor i=1:numel(filenamesAll)
    try
        [cp,ct,p]=analysisTLS_TLS(filenamesAll(i), p_Protein, maskType, cellRegion, verbose);
    catch ME
        disp(['ERROR: ' filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]);
        errFiles=[errFiles; filenamesAll(i)];
        errFileNames=[errFileNames, [filenamesAll(i).baseName num2str(filenamesAll(i).fileNumber)]];
        continue
    end
end
delete(gcp);
save(fullfile(rootFolder, 'errFiles.mat'), 'errFiles','errFileNames');


function filenamesAll=defineFileNames(rootFolder, outlineFolder, resultFolder, imageFolder, ProteinChannel, RNAChannel)
filenamesAll=[];

fileFunc=@(channel, base, num) [channel,base,num2str(num)]; %this needs to be modified according to the naming convention
baseName='ST_MS2_';
fileNumbers={2 3 4 5 6 7 8 9 10 11 12 13 14 15 16};
filenames=struct('rootFolder', rootFolder, 'outlineFolder', outlineFolder, ...
    'resultFolder', resultFolder, 'imageFolder', imageFolder, 'ProteinChannel', ProteinChannel, ...
    'RNAChannel', RNAChannel, 'baseName', baseName, 'fileFunc', fileFunc, 'fileNumber', fileNumbers);
filenamesAll=[filenamesAll, filenames];

end