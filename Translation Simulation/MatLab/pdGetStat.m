clearvars
if ispc % Running in Windows
    FDR.Base = 'E:\OneDrive - Johns Hopkins\Ha-JKwon\Projects\01 Translation Dynamics\Phase diagram\';
    FDR.Raw = 'E:\OneDrive - Johns Hopkins\Ha-JKwon\Projects\01 Translation Dynamics\Phase diagram\Raw\';
    FDR.Stat = 'E:\OneDrive - Johns Hopkins\Ha-JKwon\Projects\01 Translation Dynamics\Phase diagram\stat\';
elseif ismac % Running in Mac
    FDR.Base = '/Users/fendori/Library/CloudStorage/OneDrive-JohnsHopkins/Ha-JKwon/Projects/01 Translation Dynamics/Phase diagram/';
    FDR.Raw = '/Users/fendori/Library/CloudStorage/OneDrive-JohnsHopkins/Ha-JKwon/Projects/01 Translation Dynamics/Phase diagram/Raw/';
    FDR.Stat = '/Users/fendori/Library/CloudStorage/OneDrive-JohnsHopkins/Ha-JKwon/Projects/01 Translation Dynamics/Phase diagram/stat/';
end

%% Bursts
load([FDR.Raw,'DataB'])
Stat = zeros(size(DataB,1),4);
for i = 1:size(DataB,1)
    Stat(i,:) = getStat(DataB.MB{i,1}/6,200);
end
StatB = array2table(Stat,'VariableNames',{'MB_Mean','MB_SEM','MB_Mean_Median','MB_STD_Median'});
save([FDR.Stat,'StatB'],'StatB')
clear i ppm Stat StatB DataB

%% Dwells
load([FDR.Raw,'DataD'])
Stat = zeros(size(DataD,1),4);
for i = 1:size(DataD,1)
    Stat(i,:) = getStat(DataD.MD{i,1}/6,200);
end
StatD = array2table(Stat,'VariableNames',{'MD_Mean','MD_SEM','MD_Mean_Median','MD_STD_Median'});
save([FDR.Stat,'StatD'],'StatD')
clear i ppm Stat StatD DataD

%% #NPs
load([FDR.Raw,'DataNP'])
Stat = zeros(size(DataNP,1),4);
for i = 1:size(DataNP,1)
    Data = DataNP.NP{i,1};
    Data(Data<0.2)=[];
    Stat(i,:) = [mean(Data),std(Data)/sqrt(length(Data)),mean(Data),std(Data)/sqrt(length(Data))];
end
StatNP = array2table(Stat,'VariableNames',{'NP_Mean','NP_SEM','NP_Mean_Rep','NP_SEM_Rep'});
save([FDR.Stat,'StatNP'],'StatNP')
clear i ppm Stat StatNP DataNP

%% Functions
function Stat = getStat(Data,nBS)
%     if isempty(Data)
%         bsmed = 0;
%     elseif length(Data)==1
%         bsmed = Data;
%     elseif length(Data)>=2
%         bsmed = bootstrp(nBS,@median,Data,'Options',statset('UseParallel',false));
%     end
%     Stat = [mean(Data),std(Data)/sqrt(length(Data)),mean(bsmed),std(bsmed)];
    Stat = [mean(Data),std(Data)/sqrt(length(Data)),median(Data),1];
end
