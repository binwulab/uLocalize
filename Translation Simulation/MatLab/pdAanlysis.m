%% Clearing and setting folders
clearvars
if ispc % Running in Windows
    FDR.base = 'E:\OneDrive - Johns Hopkins\Ha-JKwon\Projects\01 Translation Dynamics\Phase diagram\';
    FDR.raw = 'E:\OneDrive - Johns Hopkins\Ha-JKwon\Projects\01 Translation Dynamics\Phase diagram\Raw\';
    FDR.stat = 'E:\OneDrive - Johns Hopkins\Ha-JKwon\Projects\01 Translation Dynamics\Phase diagram\stat\';
elseif ismac % Running in Mac
    FDR.base = '/Users/fendori/Library/CloudStorage/OneDrive-JohnsHopkins/Ha-JKwon/Projects/01 Translation Dynamics/Phase diagram/';
    FDR.raw = '/Users/fendori/Library/CloudStorage/OneDrive-JohnsHopkins/Ha-JKwon/Projects/01 Translation Dynamics/Phase diagram/Raw/';
    FDR.stat = '/Users/fendori/Library/CloudStorage/OneDrive-JohnsHopkins/Ha-JKwon/Projects/01 Translation Dynamics/Phase diagram/stat/';
end

%% Input arguments
arginMeanMed = 'Mean'; % 'Mean' or 'Median'
arginCI = 1;
arginCI_Intersect = {'Mean';0.2}; % arginCI_Intersect{1} = 'SEM' or 'Mean'

%% Script
[IN,N,PD,resExp,arrayExp,arrayL] = loadVars(arginMeanMed,FDR);
PDnRMSE = pdGetRMSE(N,PD,arrayExp,arrayL); % PDnRMSE = sqrt(meanPD^2+meanExp^2)/semExp
resFit = pdFitExp(arginCI,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp);
resGFit = pdFitGlobal(arginCI,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp);
% resGFit = pdGetCI_Boundary(arginCI,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp,resGFit);
[resFit,resGFit] = pdGetCI_Intersect(arginCI_Intersect,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp,resFit,resGFit); % Get CI from the intersects of each measurement error

%% Figures
% figure(1)
%     [x,y] = meshgrid(IN.ID,IN.IB);
%     bd = table2array(resGFit.Boundary.ORF.IBID_MB);
%     offset = 10;
%     surf(x,y,squeeze(PD(16,:,:,1,1)),'EdgeColor','none'); view(2)
%     hold on
%     plot3(x(1,28),y(13,1),(PD(16,13,28,1,1))+offset,'.k','MarkerSize',10);
%     plot3(IN.ID(bd(:,1)),IN.IB(bd(:,2)),bd(:,3)+offset,'-k','LineWidth',1);
%     hold off
%     axis tight
% clear x y bd offset

%% Functions
function [IN,N,PD,resExp,arrayExp,arrayL] = loadVars(arginMeanMed,FDR)
    load([FDR.base,'DataStat.mat'],'IN','N')
    load([FDR.base,'PD.mat'],'PD')
    if strcmp(arginMeanMed,'Mean')
        PD = PD(:,:,:,:,[1,2,5,6,9,10]);
        load([FDR.base,'ResExpMean.mat'],'resExp','arrayExp','arrayL')
    elseif strcmp(arginMeanMed,'Median')
        PD = PD(:,:,:,:,[3,4,7,8,11,12]);
        load([FDR.base,'ResExpMedian.mat'],'resExp','arrayExp','arrayL')
    end
end