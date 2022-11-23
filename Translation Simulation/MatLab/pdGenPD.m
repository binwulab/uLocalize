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

load([FDR.raw,'DataIN'],'DataIN','IN','N');
load([FDR.stat,'StatB']); StatB.Properties.VariableNames = {'MB_Mean','MB_SEM','MB_Mean_Median','MB_STD_Median'};
load([FDR.stat,'StatD']); StatD.Properties.VariableNames = {'MD_Mean','MD_SEM','MD_Mean_Median','MD_STD_Median'};
load([FDR.stat,'StatNP']); StatNP.Properties.VariableNames = {'NP_Mean','NP_SEM','NP_Mean_Rep','NP_SEM_Rep'};  

DataStat = [DataIN(:,1:4),StatB,StatD,StatNP];
clear DataIN StatB StatD StatNP

ArrayStat = table2array(DataStat);   
PD = zeros(N.ki,N.koff,N.kon,N.L,12); % 12 = { 3 measurables (MB, MD, NP) } x { 4 quantities (Mean, SEM, Mean_Median, STD_Median) }
for i = 1:N.ki
    for j = 1:N.koff
        for k = 1:N.kon
            for l = 1:N.L
               idx = find(ArrayStat(:,1)==IN.ki(i)&ArrayStat(:,2)==IN.koff(j)&ArrayStat(:,3)==IN.kon(k)&ArrayStat(:,4)==IN.L(l));
               if ~isempty(idx)
                   for m = 1:12
                       PD(i,j,k,l,m) = ArrayStat(idx,m+4);
                   end
               end
            end
        end
    end
end

save([FDR.base,'PD'],'PD')
save([FDR.base,'DataStat'],'DataStat','IN','N')