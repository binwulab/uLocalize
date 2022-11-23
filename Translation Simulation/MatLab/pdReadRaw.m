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

%% Read data sequencially
for Mode = 1:4
    tic
    %--- Read designated result
    Data01 = pdReadRawCSV([FDR.Raw,'res0-35.csv'],Mode); % 1st, all, 500 mRNA
    Data02 = pdReadRawCSV([FDR.Raw,'res36-67.csv'],Mode); % 1st, all, 500 mRNA
    Data03 = pdReadRawCSV([FDR.Raw,'res68-99.csv'],Mode); % 1st, all, 500 mRNA
    Data04 = pdReadRawCSV([FDR.Raw,'res-missing-983.csv'],Mode); % 1st, all, 500 mRNA
    Data05 = pdReadRawCSV([FDR.Raw,'res100-199.csv'],Mode); % 1st, all, 500 mRNA
    Data06 = pdReadRawCSV([FDR.Raw,'res200-299.csv'],Mode); % 1st, all, 500 mRNA
    Data07 = pdReadRawCSV([FDR.Raw,'res300-399.csv'],Mode); % 1st, all, 500 mRNA
    Data08 = pdReadRawCSV([FDR.Raw,'BW_res.csv'],Mode); % 1st, all, 500 mRNA
    Data09 = pdReadRawCSV([FDR.Raw,'BW_res_01.csv'],Mode); % 1st, all, 500 mRNA
    Data10 = pdReadRawCSV([FDR.Raw,'BW_res_02.csv'],Mode); % 1st, all, 500 mRNA
    Data11 = pdReadRawCSV([FDR.Raw,'BW_res_220710.csv'],Mode); % 1st, all, 500 mRNA
    Data12 = pdReadRawCSV([FDR.Raw,'res-080322.csv'],Mode); % 2nd, all, 500 mRNA
    Data51 = pdReadRawCSV([FDR.Raw,'r-long-burst.csv'],Mode); % 1st, IB>30, 1000 mRNA
    Data52 = pdReadRawCSV([FDR.Raw,'r-long-burst-2.csv'],Mode); % 1st, IB>30, 1000 mRNA
    Data53 = pdReadRawCSV([FDR.Raw,'res-080622.csv'],Mode); % 3rd, all, 1000 mRNA
    Data54 = pdReadRawCSV([FDR.Raw,'res-081022.csv'],Mode); % 4th, all, 1000 mRNA
    Data55 = pdReadRawCSV([FDR.Raw,'res-081922.csv'],Mode); % 5th, long IBs (5 conditions), 1000 mRNA
    Data56 = pdReadRawCSV([FDR.Raw,'res-082622-pt1.csv'],Mode); % 5th, long IBs (5 conditions), 1000 mRNA
    Data57 = pdReadRawCSV([FDR.Raw,'res-082622-pt2.csv'],Mode); % 5th, long IBs (5 conditions), 1000 mRNA
    Data58 = pdReadRawCSV([FDR.Raw,'res-082622-pt3.csv'],Mode); % 5th, long IBs (5 conditions), 1000 mRNA
    %--- Merge Data into a string array and generate nRNA information
    Raw500 = table2array([Data01;Data02;Data03;Data04;Data05;Data06;Data07;Data08;Data09;Data10;Data11;Data12]);
    nRNA500 = 500*ones(size(Raw500,1),1);
    Raw1000 = table2array([Data51;Data52;Data53;Data54;Data55;Data56;Data57;Data58]);
    nRNA1000 = 1000*ones(size(Raw1000,1),1);
    clear Data01 Data02 Data03 Data04 Data05 Data06 Data07 Data08 Data09 Data10 Data11 Data12
    clear Data51 Data52 Data53 Data54 Data55 Data56 Data57 Data58
    RawINBDNP = [Raw500;Raw1000];
    nRNA = [nRNA500;nRNA1000];
    clear Raw500 nRNA500 Raw1000 nRNA1000
    if Mode == 1
        %--- Find and remove dummy Data
        IdxDummy = find(RawINBDNP=="kin,koff,kon,len");
        RawINBDNP(IdxDummy) = [];
        nRNA(IdxDummy) = [];
        %--- Generate a numeric array of inputs
        NumIN = zeros(size(RawINBDNP,1),5);
        for i = 1:size(RawINBDNP,1)
            var = str2num(RawINBDNP(i));
            NumIN(i,1:4) = var(1:4);
        end
        NumIN(:,5) = nRNA;
        clear i var RawINBDNP nRNA
        %--- Merge very close inputs // only for dataset 1
        NumIN(NumIN(:,3)==0.009980,3) = 0.010000;
        NumIN(NumIN(:,3)==0.012255,3) = 0.012222;
        NumIN(NumIN(:,2)==0.001400,2) = 0.001401;
        %--- Get input parameter information
        IN.ki = unique(NumIN(:,1)); N.ki = size(IN.ki,1);
        IN.koff = unique(NumIN(:,2)); N.koff = size(IN.koff,1);
        IN.kon = unique(NumIN(:,3)); N.kon = size(IN.kon,1);
        IN.L = unique(NumIN(:,4)); N.L = size(IN.L,1);
        %--- Remove undesired inputs // only for dataset 1
        IN.ki(28:end) = []; N.ki = size(IN.ki,1);
        IN.koff([1,3,5,20,23,26,28,30,32,38,40,41,43:45]) = []; N.koff = size(IN.koff,1);
        IN.kon([33,35,37,39:41,43:end]) = []; N.kon = size(IN.kon,1);
        IdxCensored = ismember(NumIN(:,1),IN.ki) & ismember(NumIN(:,2),IN.koff) & ismember(NumIN(:,3),IN.kon);
        NumIN = NumIN(IdxCensored,:);
        %--- Add input information in alternative units
        IN.nInit = IN.ki*60; IN.IB = 1./IN.koff/60; IN.ID = 1./IN.kon/60;
        %--- Merge replicates
        var = strings(size(NumIN,1),1);
        for i = 1:size(NumIN,1)
            var(i) = strcat(num2str(NumIN(i,1)),',',num2str(NumIN(i,2)),',',num2str(NumIN(i,3)),',',num2str(NumIN(i,4)));
        end
        [~,IdxFirst] = unique(var,'first');
        IdxDup = find(not(ismember((1:numel(var))',IdxFirst)));
        IdxTarget = nan(size(IdxDup,1),1);
        for i = 1:size(IdxDup,1)
           Idx = find(var==var(IdxDup(i)));
           IdxTarget(i) = Idx(1);
           NumIN(IdxTarget(i),5) = NumIN(IdxTarget(i),5)+NumIN(IdxDup(i),5);
           NumIN(IdxDup(i),:) = nan;
        end
        NumIN(isnan(NumIN(:,1)),:) = [];
        clear var i Idx IdxFirst
        %--- Convert result into a table
        DataIN = array2table(NumIN,'VariableNames',{'ki','koff','kon','L','nRNA'});
        save([FDR.Raw,'DataIN'],'DataIN','IdxCensored','IdxDummy','IdxDup','IdxTarget','IN','N')
        clear NumIN DataIN
    else
        load([FDR.Raw,'DataIN'],'IdxCensored','IdxDummy','IdxDup','IdxTarget','IN','N')
        %--- Remove dummy
        RawINBDNP(IdxDummy) = [];
        %--- Remove undesired bursts
        RawINBDNP = RawINBDNP(IdxCensored,:);
        %--- Merge replicates
        for i = 1:size(IdxDup,1)
           RawINBDNP(IdxTarget(i)) = strcat(RawINBDNP(IdxTarget(i)),',',RawINBDNP(IdxDup(i)));
           RawINBDNP(IdxDup(i)) = missing;
        end
        RawINBDNP(ismissing(RawINBDNP)) = [];
        clear i
        %--- Convert string array to a table with numeric values
        NumBDNP = cell(size(RawINBDNP,1),1);
        for i = 1:size(RawINBDNP,1)
            NumBDNP{i,1} = str2num(RawINBDNP(i));
        end
        if Mode == 2
            DataB = cell2table(NumBDNP,'VariableNames',{'MB'});
            save([FDR.Raw,'DataB'],'DataB')
            clear i RawINBDNP NumBDNP DataB
        elseif Mode == 3
            DataD = cell2table(NumBDNP,'VariableNames',{'MD'});
            save([FDR.Raw,'DataD'],'DataD')
            clear i RawINBDNP NumBDNP DataD
        elseif Mode == 4
            DataNP = cell2table(NumBDNP,'VariableNames',{'NP'});
            save([FDR.Raw,'DataNP'],'DataNP')
            clear i RawINBDNP NumBDNP DataNP
        end
    end
    toc
end

