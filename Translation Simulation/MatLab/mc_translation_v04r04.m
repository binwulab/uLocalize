%% Notes
%--- v04r04: Updated on June 16, 2022
% Implementing single-line parfor loop
%--- v04r03: Updated on June 08, 2022
% Implementing cell-based analysis
% Introducing Param.Frame.CO, fCO
%--- v04r02: Updated on May 26, 2022
% Implementing Gamma-distributed burst widths
%--- v04r01: Updated on MAY 04, 2022; |CAUTION|Not compatible with all previous versions
% Re-design all the variables to avoid 'structure' or 'cell' during the intensive loop-based calculation

%% Mode selection
Mode.MC = 1; % Monte Carlo
if Mode.MC == 1; clearvars; Mode.MC = 1;
    Mode.E = 3; % Elongation; |0|SeqDet |1|RndDet |2|SeqTASEP |3|RndTASEP
    Mode.B = 1; % Bursting; |0|Constitutive |1|2SRTwOvlp |2|2SRTwoOvlp
    Mode.LT = 1; % Track length (LT); |0|FixLT |1|ExpLT
    Mode.FRAP = 0; % FRAP; |0|NoFRAP |1|FRAP
    Mode.RO = 0; % Run-off; |0|NoRO |1|RO    
end
Mode.nR0 = 1; % Population of zero ribosome; |0|Preserve |1|Remove
Mode.nRC = 0; % Correction of #R; |0|Raw |1|Corrected
Mode.Fig = 1; % Generating figures; |0|No operation |1|Generating figures
Mode.ISP = 0; % Inspection of FI traces; |0|No action |1|Check FI traces molecule-by-molecule

%% Parameters
if Mode.MC == 1
    %--- #RNA
    N = 4000;
    Param.Cell.NpC = 100; % #RNA per cell
    %--- RNA reporter information; AUG-ORF1-Tag1-ORF2-polyA-Tag2-ORF3-STOP
    Param.RNA.l1 = 0; % Length of ORF1 (aa)
    Param.RNA.nt1 = 24; % #Tag1
    Param.RNA.lt1 = 24; % Length of each Tag1 (aa)
    Param.RNA.l2 = 265; % Length of ORF2 (aa); |265|841 |506|1082 |678|1254 |817|1393 |1063|1639
    Param.RNA.lpA = 0; % Length of polyA (aa)
    Param.RNA.nt2 = 0; % #Tag2
    Param.RNA.lt2 = 0; % Length of each Tag2 (aa)
    Param.RNA.l3 = 0; % Length of ORF3                                  
    %--- Kinetic rates
    Param.Rate.i = 0.0016667;%0.080603887234854; % Initiation rate (s-1)
    Param.Rate.e = 4.7; % Elongation rate (aa/s)
    Param.Rate.eA = 0.5*Param.Rate.e; % Elongation rate on polyA (aa/s)
    Param.Rate.t = 0.5*Param.Rate.e; % Termination rate (s-1)
    Param.Rate.on = 0.005128; % On-switching rate (s-1)
    Param.Rate.off = 0.000347; % Off-switching rate (s-1)
    %--- Timing parameters
    Param.Time.Exp = 10; % Experimental exposure time (unit of frame) (s)
    if Mode.E <= 1; Param.Time.dT = 0.2; else; Param.Time.dT = 0.02; end % Simulation dt; determined by the elongation mode
    %--- Frame parameters
    Param.Frame.LTmin = 180; % Minimum LT
    Param.Frame.LTmu = 355.12; % Mu of Exp-distributed LT
    Param.Frame.LTmax = 1081; % Maximum LT
    Param.Frame.Eq = 100; % Buffer frame for equilibrium
    Param.Frame.FRAP = 60; % Frame for FRAP
    Param.Frame.RO = 60; % Frame for RO
    Param.Frame.CO = 3; % Cut-off frame during binarization
    %--- Photon parameters; Gamma-distributed photon numbers for single tag
    Param.Photon.shape = 0.27;
    Param.Photon.scale = 236.2;
    Param.Photon.tag = Param.Photon.shape*Param.Photon.scale;
    %--- Ribosome footprint
    Param.RFP = 9; % In 'aa'
    %--- Shape of gamma-distributed burst and initiation times
    Param.BTS = 1; % |1|Exponential |N|Gamma
    Param.ITS = 1; % |1|Exponential |N|Gamma
    %--- Calculated parameters   
    L = Param.RNA.l1+Param.RNA.nt1*Param.RNA.lt1+Param.RNA.l2+Param.RNA.lpA+Param.RNA.nt2*Param.RNA.lt2+Param.RNA.l3; % Length of RNA
end

%% Decompose variables for parallel computation
modeE = Mode.E; modeB = Mode.B; modeLT = Mode.LT; modeFRAP = Mode.FRAP; modeRO = Mode.RO; modenR0 = Mode.nR0; modenRC = Mode.nRC;
ki = Param.Rate.i; ke = Param.Rate.e; keA = Param.Rate.eA; kt = Param.Rate.t; kon = Param.Rate.on; koff = Param.Rate.off;
LTmin = Param.Frame.LTmin; LTmu = Param.Frame.LTmu; LTmax = Param.Frame.LTmax; fEq = Param.Frame.Eq; fFRAP = Param.Frame.FRAP; fRO = Param.Frame.RO; fCO = Param.Frame.CO;
tdT = Param.Time.dT; tExp = Param.Time.Exp;
l1 = Param.RNA.l1; nt1 = Param.RNA.nt1; lt1 = Param.RNA.lt1; l2 = Param.RNA.l2; lpA = Param.RNA.lpA; nt2 = Param.RNA.nt2; lt2 = Param.RNA.lt2; l3 = Param.RNA.l3;
pShape = Param.Photon.shape; pScale = Param.Photon.scale; pTag = Param.Photon.tag;
RFP = Param.RFP; BTS = Param.BTS; ITS = Param.ITS;  

%% Monte Carlo
if Mode.MC == 1
    Data(N) = struct('LT',[],'nS',[],'TSI',[],'ICD',[],'nR',[],'nT1',[],'nT2',[],'FIT1',[],'FIT1B',[],...
        'FIT2',[],'CD',[],'CO',[],'iR',[],'iMP',[],'nRS',[],'FIT1S',[],'BI',[],'BM',[]);
    %--- Activate parallel pool and set progress bar
    getPool = gcp('nocreate');
    if isempty(getPool); parpool(feature('NumThreads')); end
    ppm = ParforProgressStarter2(' Simulating... ',N,0.1,1,1,1);
    %--- Simulation for individual RNAs; 'i' = RNA index
    parfor i = 1:N        
        Data(i) = runmc(Mode,Param,L);
        %--- Update progress bar
        ppm.increment(i);
    end
    %--- Clear variables
    delete(ppm);
    clear ppm getPool
end

%% Analysis
%--- Initialize variables
nRraw = []; nRT1 = [];
IT = []; TT = [];
CO = [];
FIFRAP = []; FIRO = [];
BI = []; DI = []; BM = []; DM = [];
%--- Individual RNA analysis
for i = 1:N
    nRS = Data(i).nRS; iMP = Data(i).iMP; FIT1 = Data(i).FIT1;
    bdi = Data(i).BI; bdm = Data(i).BM;
    %--- #Ribosome
    nRFIT1S = Data(i).FIT1S/(nt1*pTag);
    if modenRC == 1; nRFIT1S = nRFIT1S*((L*nt1)/((L-nt1*lt1)*nt1+lt1*sum(0:1:nt1-1))); end
%     nRFIT1S = round(nRFIT1S);
    if modenR0 == 1; nRS(nRS==0) = []; nRFIT1S(nRFIT1S<=0.2) = []; end
    nRraw = [nRraw;[i*ones(length(nRS),1),nRS]]; nRT1 = [nRT1;[i*ones(length(nRFIT1S),1),nRFIT1S]];
    %--- Initiation time
    if ~isempty(iMP)
    sIT = iMP(:,2); sIT(2:end) = sIT(2:end)-sIT(1:end-1); sIT(1) = []; IT = [IT;sIT*tdT];
    end
    %--- Translation time
    if ~isempty(iMP)
    sTT = iMP(:,11)-iMP(:,2); TT = [TT;sTT*tdT]; 
    end
    %--- Codon occupancy
    CO(i,:) = Data(i).CO;
    %--- FRAP
    if modeFRAP == 1; FIFRAP(i,:) = FIT1(end-fFRAP:end)/mean(FIT1(1:end-(fFRAP+1))); end
    %--- RO
    if modeRO == 1; FIRO(i,:) = FIT1(end-fRO:end)/mean(FIT1(1:end-(fRO+1))); end
    %--- Bursts and dwells
    idx = find(bdi(:,3)==1); BI = [BI;[i*ones(length(idx),1),bdi(idx,4)]];
    idx = find(bdi(:,3)==0); DI = [DI;[i*ones(length(idx),1),bdi(idx,4)]];
    idx = find(bdm(:,3)==1); BM = [BM;[i*ones(length(idx),1),bdm(idx,4)]];
    idx = find(bdm(:,3)==0); DM = [DM;[i*ones(length(idx),1),bdm(idx,4)]];
end
if N > 1
    idx = isnan(CO(:,1)); CO(idx==1,:) = []; CO = mean(CO);
    if modeFRAP == 1; FRAP.FI = FIFRAP; FRAP.FI_Mean = mean(FIFRAP); FRAP.FI_Se = std(FIFRAP)/sqrt(N); end
    if modeRO == 1; RO.FI = FIRO; RO.FI_Mean = mean(FIRO); RO.FI_Se = std(FIRO)/sqrt(N); end
end
CD = [Data.CD]';
if modeFRAP == 1
    FRAP.Time = (0:1:fFRAP)*tExp;
    FRAP.Theory = ((ke/nt1)*FRAP.Time)/(round(L/lt1)-((nt1-1)/2));
    FRAP.Theory(FRAP.Theory>1) = 1;
end
if modeRO == 1
    RO.Time = (0:1:fRO)*tExp;
    RO.Theory = f_RO_FI([ke;0],RO.Time',[nt1;lt1;l2])';
end
%--- Summarized into structures
Ribosome.nRraw = nRraw; Ribosome.nRT1 = nRT1; Ribosome.IT = IT; Ribosome.TT = TT; Ribosome.CO_Mean = CO; Ribosome.CD = CD;
BD.BI = BI; BD.DI = DI; BD.BM = BM; BD.DM = DM;
%--- Get some values
Value.nR_Mean_Theory = ki*L/ke; Value.nR_Mean = mean(nRraw(:,2)); Value.nRT1_Mean = mean(nRT1(:,2));
Value.IT_Mean = mean(IT); Value.TT_Mean = mean(TT);
Value.FI_Mean_Theory =nt1*pTag*(ki/(ke/lt1))*(L/nt1-(1/2)*(nt1));
Value.BI_Mean = mean(BD.BI(:,2)); Value.BI_Median = median(BD.BI(:,2)); Value.DI_Mean = mean(BD.DI(:,2)); Value.DI_Median = median(BD.DI(:,2));
Value.BM_Mean = mean(BD.BM(:,2)); Value.BM_Median = median(BD.BM(:,2)); Value.DM_Mean = mean(BD.DM(:,2)); Value.DM_Median = median(BD.DM(:,2));
Value.LT_Mean = mean([Data.LT]); Value.LT_Median = median([Data.LT]);
%--- Cell-based analysis

Param.Cell.NC = ceil(N/Param.Cell.NpC);
Cell = struct();
NC = Param.Cell.NC;
NpC = Param.Cell.NpC; 
for i = 1:NC
    %--- Get information per cell   
    lb = (i-1)*NpC+1; ub = i*NpC; if i == NC; ub = N; end
    Cell(i).nRraw = nRraw(nRraw(:,1)>=lb&nRraw(:,1)<=ub,2);
    Cell(i).nRrawMean = mean(nRraw(nRraw(:,1)>=lb&nRraw(:,1)<=ub,2));
    Cell(i).nRT1 = nRT1(nRT1(:,1)>=lb&nRT1(:,1)<=ub,2);
    Cell(i).nRT1Mean = mean(nRT1(nRT1(:,1)>=lb&nRT1(:,1)<=ub,2));
    Cell(i).BI = BI(BI(:,1)>=lb&BI(:,1)<=ub,2);
    Cell(i).BIMean = mean(BI(BI(:,1)>=lb&BI(:,1)<=ub,2));
    Cell(i).BIMedian = median(BI(BI(:,1)>=lb&BI(:,1)<=ub,2));
    Cell(i).DI = DI(DI(:,1)>=lb&DI(:,1)<=ub,2);
    Cell(i).DIMean = mean(DI(DI(:,1)>=lb&DI(:,1)<=ub,2));
    Cell(i).DIMedian = median(DI(DI(:,1)>=lb&DI(:,1)<=ub,2));
    Cell(i).BM = BM(BM(:,1)>=lb&BM(:,1)<=ub,2);
    Cell(i).BMMean = mean(BM(BM(:,1)>=lb&BM(:,1)<=ub,2));
    Cell(i).BMMedian = median(BM(BM(:,1)>=lb&BM(:,1)<=ub,2));
    Cell(i).DM = DM(DM(:,1)>=lb&DM(:,1)<=ub,2);
    Cell(i).DMMean = mean(DM(DM(:,1)>=lb&DM(:,1)<=ub,2));
    Cell(i).DMMedian = median(DM(DM(:,1)>=lb&DM(:,1)<=ub,2));
    Cell(i).LT = [Data(lb:ub).LT]';
    if modeB > 0
        %--- Generate and fit CDFs of burst, dwell, and LT per cell
        figure(10)
        fitopt = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'Display','off');
        X = [1:1:LTmax]-0.5; Cell(i).CDFBin = X;
        temphist = histogram(BM(BM(:,1)>=lb&BM(:,1)<=ub,2),'Normalization','CDF','BinWidth',1,'BinLimits',[0,LTmax],'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
        Y = temphist.Values; Cell(i).CDFBM = Y; iCDFBM = 1-Y;
        fitfunc = @(x)((1-exp(-X/x))-Y);
        [res,~] = lsqnonlin(fitfunc,Value.BM_Mean,[],[],fitopt); Cell(i).BMFitMean = res;
        temphist = histogram(DM(DM(:,1)>=lb&DM(:,1)<=ub,2),'Normalization','CDF','BinWidth',1,'BinLimits',[0,LTmax],'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
        Y = temphist.Values; Cell(i).CDFDM = Y; iCDFDM = 1-Y;
        fitfunc = @(x)((1-exp(-X/x))-Y);
        [res,~] = lsqnonlin(fitfunc,Value.DM_Mean,[],[],fitopt); Cell(i).DMFitMean = res;
        temphist = histogram([Data(lb:ub).LT],'Normalization','CDF','BinWidth',1,'BinLimits',[0,LTmax],'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
        Y = temphist.Values; Cell(i).CDFLT = Y;
        fitfunc = @(x)(f_LT_CDF(X',LTmin,x,LTmax)-Y');
        [res,~] = lsqnonlin(fitfunc,LTmu,[],[],fitopt); Cell(i).LTFitMean = res; iCDFLT = 1-expcdf(X,LTmin+res);
        %--- LT correction and deduce intrinsic values
        iY = iCDFBM./iCDFLT; Y = 1-iY; Cell(i).CDFBMCor = Y;
        fitfunc = @(x)((1-exp(-X/x))-Y);
        [BMCorFitMean,~] = lsqnonlin(fitfunc,mean(Y),[],[],fitopt); Cell(i).BMCorFitMean = BMCorFitMean;
        iY = iCDFDM./iCDFLT; Y = 1-iY; Cell(i).CDFDMCor = Y;
        fitfunc = @(x)((1-exp(-X/x))-Y);
        [DMCorFitMean,~] = lsqnonlin(fitfunc,mean(Y),[],[],fitopt); Cell(i).DMCorFitMean = DMCorFitMean;
        konCal = 1/((DMCorFitMean-fCO)*tExp); Cell(i).konCal = konCal; POV = 1-exp(-konCal*L/ke);
        koffCal = konCal/(BMCorFitMean*tExp*konCal*(1-POV)+POV*konCal*(L/ke)-POV-kon*(L/ke)); Cell(i).koffCal = koffCal;
        close figure 10
    end
end
%--- Summarize values
CellValue.nR = [mean([Cell.nRrawMean]),std([Cell.nRrawMean])]; CellValue.nRT1 = [mean([Cell.nRT1Mean]),std([Cell.nRT1Mean])];
CellValue.BI = [mean([Cell.BIMean]),std([Cell.BIMean]),mean([Cell.BIMedian]),std([Cell.BIMedian])]; 
CellValue.DI = [mean([Cell.DIMean]),std([Cell.DIMean]),mean([Cell.DIMedian]),std([Cell.DIMedian])]; 
CellValue.BM = [mean([Cell.BMMean]),std([Cell.BMMean]),mean([Cell.BMMedian]),std([Cell.BMMedian])];
CellValue.DM = [mean([Cell.DMMean]),std([Cell.DMMean]),mean([Cell.DMMedian]),std([Cell.DMMedian])];
if modeB > 0
    CellValue.BM = [mean([Cell.BMMean]),std([Cell.BMMean]),mean([Cell.BMMedian]),std([Cell.BMMedian]),mean([Cell.BMFitMean]),std([Cell.BMFitMean]),mean([Cell.BMCorFitMean]),std([Cell.BMCorFitMean])];
    CellValue.DM = [mean([Cell.DMMean]),std([Cell.DMMean]),mean([Cell.DMMedian]),std([Cell.DMMedian]),mean([Cell.DMFitMean]),std([Cell.DMFitMean]),mean([Cell.DMCorFitMean]),std([Cell.DMCorFitMean])];
    CellValue.CDFBin = Cell(1).CDFBin';
    CellValue.CDFBM = [mean(reshape([Cell.CDFBM],[LTmax,NC]),2),std(reshape([Cell.CDFBM],[LTmax,NC]),1,2)/sqrt(NC)];
    CellValue.CDFDM = [mean(reshape([Cell.CDFDM],[LTmax,NC]),2),std(reshape([Cell.CDFDM],[LTmax,NC]),1,2)/sqrt(NC)];
    CellValue.CDFLT = [mean(reshape([Cell.CDFLT],[LTmax,NC]),2),std(reshape([Cell.CDFLT],[LTmax,NC]),1,2)/sqrt(NC)];
    CellValue.CDFBMCor = [mean(reshape([Cell.CDFBMCor],[LTmax,NC]),2),std(reshape([Cell.CDFBMCor],[LTmax,NC]),1,2)/sqrt(NC)];
    CellValue.CDFDMCor = [mean(reshape([Cell.CDFDMCor],[LTmax,NC]),2),std(reshape([Cell.CDFDMCor],[LTmax,NC]),1,2)/sqrt(NC)];
    CellValue.konCal = [mean([Cell.konCal]),std([Cell.konCal])]; CellValue.koffCal = [mean([Cell.koffCal]),std([Cell.koffCal])];
end

%% Clean up temporary variables
clear modeE modeB modeLT modeFRAP modeRO modenR0 modenRC
clear ki ke keA kt kon koff c_ke
clear LTmin LTmu LTmax fEq fFRAP fRO fCO
clear tdT tExp
clear l1 nt1 lt1 l2 lpA nt2 lt2 l3 c_nT1 c_nT2
clear pShape pScale pTag
clear RFP BTS ITS getPool
clear nRraw nRT1 nRS nRFIT1S sIT sTT
clear IT TT CO CD FIT1 FIFRAP FIRO BI DI BM DM NC NpC POV Y
clear BMCorFitMean DMCorFitMean fitfunc fitopt iCDFBM iCDFDM iCDFLT iMP iY
clear idx bdi bdm koffCal konCal lb res temphist ub

%% Graphics
if Mode.Fig == 1
    i = 1;
    Fig.n = figure(1);
    Fig.MP = get(0, 'MonitorPositions');
    Fig.n.Position = Fig.MP(1,4)*[0.1,0.1,0.8,0.7];
    Fig.group = uitabgroup(Fig.n);
    Fig.tab1 = uitab(Fig.group, 'Title', 'TimeTrace');
    Fig.axes1 = axes('Parent', Fig.tab1);
        subplot(3,1,1)
            plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).nR,'Color','#666666','LineWidth',1)
            hold on; plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).TSI*max(Data(i).nR),'-k','LineWidth',1); hold off
            title('ITS (black) and #R (gray)'); ylabel('#R'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(Data(i).nR))
        subplot(3,1,2)
            plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).nT1,'Color','#666666','LineWidth',1)
            hold on; plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).nT2,'Color','#999999','LineWidth',1); hold off
            title('#T1 (gray) and #T2 (light gray)'); ylabel('#Tags'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(Data(i).nT1))
        subplot(3,1,3)
            plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).FIT1,'Color','#666666','LineWidth',1)
            hold on
            plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).FIT2,'Color','#999999','LineWidth',1)
            plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).FIT1B*max(Data(i).FIT1),'Color','#000000','LineWidth',1)
            if Mode.FRAP == 1
                plot([1,Data(i).LT-Param.Frame.FRAP-1]/(60/Param.Time.Exp),[1,1]*Value.FI_Mean_Theory,':r','LineWidth',1);
                plot((Data(i).LT-Param.Frame.FRAP:1:Data(i).LT)/(60/Param.Time.Exp),FRAP.Theory*Value.FI_Mean_Theory,'--r','LineWidth',1);
            end
            if Mode.RO == 1
                plot([1,Data(i).LT-Param.Frame.RO-1]/(60/Param.Time.Exp),[1,1]*Value.FI_Mean_Theory,':r','LineWidth',1);
                plot((Data(i).LT-Param.Frame.RO:1:Data(i).LT)/(60/Param.Time.Exp),RO.Theory*Value.FI_Mean_Theory,'--r','LineWidth',1);
            end
            hold off
            title('FIT1 (gray), FIT2 (light gray), and MTS (black)'); xlabel('Time (min)'); ylabel('FI (AU)'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(Data(i).FIT1))
    Fig.tab2 = uitab(Fig.group, 'Title', 'Ribosome');
    Fig.axes2 = axes('Parent', Fig.tab2);
        subplot(2,2,1)
            Hist.nR = histogram(Ribosome.nRraw(:,2),'Normalization','Probability','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            hold on
            Hist.nRT1 = histogram(Ribosome.nRT1(:,2),'Normalization','Probability','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#000000','LineWidth',1);
            hold off
            legend(num2str(Value.nR_Mean,'%.1f'),num2str(Value.nRT1_Mean,'%.1f'))
            title('#R (gray) and #R-T1 (black)'); xlabel('#Ribosome'); ylabel('Frequency'); ylim([-0.02,1.2]*max(max(Hist.nR.Values),max(Hist.nRT1.Values)))
        subplot(2,2,2)
            Hist.IT = histogram(Ribosome.IT/60,'Normalization','Probability','BinWidth',0.05,'BinLimits',[0,ceil(L/Param.Rate.e/60)],'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            hold on
            Hist.TT = histogram(Ribosome.TT/60,'Normalization','Probability','BinWidth',0.05,'DisplayStyle','Stairs','EdgeColor','#000000','LineWidth',1);
            hold off
            legend([num2str(Value.IT_Mean,'%.1f'),' s'],[num2str(Value.TT_Mean,'%.1f'),' s'])
            title('Initiation (gray) & Translation (black) times'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*ceil(Hist.TT.BinLimits(2))); ylim([-0.02,1.2]*max(max(Hist.IT.Values),max(Hist.TT.Values)))
        subplot(2,2,3)
            plot((1:L),Ribosome.CO_Mean,'Color','#666666','LineWidth',1)
            title('Codon occupancy'); xlabel('Codon position'); ylabel('Occupancy'); xlim([-0.02,1.02]*L); ylim([-0.02,1.2]*max(Ribosome.CO_Mean))
        subplot(2,2,4)
            Hist.CD = histogram(Ribosome.CD,'Normalization','Probability','DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            title('Collision density'); xlabel('Density (/frame/codon)'); ylabel('Frequency'); ylim([-0.02,1.2]*max(Hist.CD.Values))
    if Mode.FRAP == 1    
        Fig.tab3 = uitab(Fig.group, 'Title', 'FRAP');
        Fig.axes3 = axes('Parent', Fig.tab3);
            if N > 1
                Fig.FRAP_Err_x = [FRAP.Time/60,fliplr(FRAP.Time/60)];
                Fig.FRAP_Err_a = fill(Fig.FRAP_Err_x,[FRAP.FI_Mean+FRAP.FI_Se,fliplr(FRAP.FI_Mean-FRAP.FI_Se)],[0.5,0.5,0.5]);
                set(Fig.FRAP_Err_a,'EdgeColor','None'); set(Fig.FRAP_Err_a,'FaceAlpha',0.5);
                hold on
                plot(FRAP.Time/60,FRAP.FI_Mean,'-k','LineWidth',1)
                plot(FRAP.Time/60,FRAP.Theory,'--r','LineWidth',1)
                hold off
            else
                plot(FRAP.Time/60,FRAP.FI_Mean,'-k','LineWidth',1)
                hold on
                plot(FRAP.Time/60,FRAP.Theory,'--r','LineWidth',1)
                hold off
            end
            title('Normalized FIT1 (black) and Theoretical FRAP (red)'); xlabel('Time (min)'); ylabel('Normalized FI (AU)'); xlim([-0.02,1.02]*Param.Frame.FRAP/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(FRAP.FI_Mean+FRAP.FI_Se))
    end
    if Mode.RO == 1    
        Fig.tab4 = uitab(Fig.group, 'Title', 'Run-off');
        Fig.axes4 = axes('Parent', Fig.tab4);
            if N > 1
                Fig.RO_Err_x = [RO.Time/60,fliplr(RO.Time/60)];
                Fig.RO_Err_a = fill(Fig.RO_Err_x,[RO.FI_Mean+RO.FI_Se,fliplr(RO.FI_Mean-RO.FI_Se)],[0.5,0.5,0.5]);
                set(Fig.RO_Err_a,'EdgeColor','None'); set(Fig.RO_Err_a,'FaceAlpha',0.5);
                hold on
                plot(RO.Time/60,RO.FI_Mean,'-k','LineWidth',1)
                plot(RO.Time/60,RO.Theory,'--r','LineWidth',1)
                hold off
            else
                plot(RO.Time/60,RO.FI_Mean,'-k','LineWidth',1)
                hold on
                plot(RO.Time/60,RO.Theory,'--r','LineWidth',1)
                hold off
            end
            title('Normalized FIT1 (black) and Theoretical run-off (red)'); xlabel('Time (min)'); ylabel('Normalized FI (AU)'); xlim([-0.02,1.02]*Param.Frame.RO/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(RO.FI_Mean+RO.FI_Se))
    end
    Fig.tab5 = uitab(Fig.group, 'Title', 'Intrinsic');
    Fig.axes5 = axes('Parent', Fig.tab5);
        subplot(2,2,1)
            Hist.BI_PDF = histogram(BD.BI(:,2)/(60/Param.Time.Exp),'Normalization','Probability','BinWidth',5,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Mean: ',num2str(Value.BI_Mean/(60/Param.Time.Exp),'%.1f'),' min'])
            title('PDF: intrinsic burst'); xlabel('Burst width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60); if ~isempty(BD.BI); ylim([-0.02,1.2]*max(Hist.BI_PDF.Values)); end
        subplot(2,2,2)
            Hist.DI_PDF = histogram(BD.DI(:,2)/(60/Param.Time.Exp),'Normalization','Probability','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Mean: ',num2str(Value.DI_Mean/(60/Param.Time.Exp),'%.1f'),' min'])
            title('PDF: intrinsic dwell'); xlabel('Dwell width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60/6); if ~isempty(BD.DI); ylim([-0.02,1.2]*max(Hist.DI_PDF.Values)); end
        subplot(2,2,3)
            Hist.BI_CDF = histogram(BD.BI(:,2)/(60/Param.Time.Exp),'Normalization','CDF','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Median: ',num2str(Value.BI_Median/(60/Param.Time.Exp),'%.1f'),' min'])
            title('CDF: intrinsic burst'); xlabel('Burst width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60); if ~isempty(BD.BI); ylim([-0.02,1.2]*max(Hist.BI_CDF.Values)); end
        subplot(2,2,4)
            Hist.DI_CDF = histogram(BD.DI(:,2)/(60/Param.Time.Exp),'Normalization','CDF','BinWidth',0.2,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Median: ',num2str(Value.DI_Median/(60/Param.Time.Exp),'%.1f'),' min'])
            title('CDF: intrinsic dwell'); xlabel('Dwell width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60/6); if ~isempty(BD.DI); ylim([-0.02,1.2]*max(Hist.DI_CDF.Values)); end
    Fig.tab6 = uitab(Fig.group, 'Title', 'Measured');
    Fig.axes6 = axes('Parent', Fig.tab6);
        subplot(2,2,1)
            Hist.BM_PDF = histogram(BD.BM(:,2)/(60/Param.Time.Exp),'Normalization','Probability','BinWidth',5,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Mean: ',num2str(Value.BM_Mean/(60/Param.Time.Exp),'%.1f'),' min'])
            title('PDF: measured burst'); xlabel('Burst width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60); if ~isempty(BD.BM); ylim([-0.02,1.2]*max(Hist.BM_PDF.Values)); end
        subplot(2,2,2)
            Hist.DM_PDF = histogram(BD.DM(:,2)/(60/Param.Time.Exp),'Normalization','Probability','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Mean: ',num2str(Value.DM_Mean/(60/Param.Time.Exp),'%.1f'),' min'])
            title('PDF: measured dwell'); xlabel('Dwell width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60/6); if ~isempty(BD.DM); ylim([-0.02,1.2]*max(Hist.DM_PDF.Values)); end
        subplot(2,2,3)
            Hist.BM_CDF = histogram(BD.BM(:,2)/(60/Param.Time.Exp),'Normalization','CDF','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Median: ',num2str(Value.BM_Median/(60/Param.Time.Exp),'%.1f'),' min'])
            title('CDF: measured burst'); xlabel('Burst width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60); if ~isempty(BD.BM); ylim([-0.02,1.2]*max(Hist.BM_CDF.Values)); end
        subplot(2,2,4)
            Hist.DM_CDF = histogram(BD.DM(:,2)/(60/Param.Time.Exp),'Normalization','CDF','BinWidth',0.2,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            legend(['Median: ',num2str(Value.DM_Median/(60/Param.Time.Exp),'%.1f'),' min'])
            title('CDF: measured dwell'); xlabel('Dwell width (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax*Param.Time.Exp/60/6); if ~isempty(BD.DM); ylim([-0.02,1.2]*max(Hist.DM_CDF.Values)); end  
    Fig.tab7 = uitab(Fig.group, 'Title', 'Track');
    Fig.axes7 = axes('Parent', Fig.tab7);
        subplot(2,2,1)
            Hist.LT_CDF = histogram([Data.LT]/(60/Param.Time.Exp),'Normalization','CDF','BinWidth',1,'DisplayStyle','Stairs','EdgeColor','#666666','LineWidth',1);
            title('CDF: track length'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(Hist.LT_CDF.Values))
    if Mode.B > 0
        Fig.tab8 = uitab(Fig.group, 'Title', 'SingleCell_1');
        Fig.axes8 = axes('Parent', Fig.tab8);
            subplot(2,2,1)
                data = [CellValue.nR(1),CellValue.nRT1(1)]; err = [CellValue.nR(2),CellValue.nRT1(2)];
                bar(data)
                hold on
                errorbar(data,err,'Color',[0,0,0],'LineStyle','None')
                dmy1 = line(nan,nan,'Linestyle','none','Marker','none','Color','none');
                dmy2 = line(nan,nan,'Linestyle','none','Marker','none','Color','none');
                hold off
                set(gca,'xticklabel',{'#R','#R-T1'});
                legend([dmy1,dmy2],{['[#R] ',num2str(CellValue.nR(1),'%.2f'),'±',num2str(CellValue.nR(2),'%.2f')],['[#R-T1] ',num2str(CellValue.nRT1(1),'%.2f'),'±',num2str(CellValue.nRT1(2),'%.2f')]}); legend('boxoff')
                title('#R and #R-T1'); ylabel('#Ribosomes'); ylim([-0.02,1.2]*max(sum(CellValue.nR),sum(CellValue.nRT1)))
            subplot(2,2,2)
                Fig.Cell_CDFLT_Err_x = [CellValue.CDFBin'/(60/Param.Time.Exp),fliplr(CellValue.CDFBin'/(60/Param.Time.Exp))];
                Fig.Cell_CDFLT_Err_a = fill(Fig.Cell_CDFLT_Err_x,[CellValue.CDFLT(:,1)'+CellValue.CDFLT(:,2)',fliplr(CellValue.CDFLT(:,1)'-CellValue.CDFLT(:,2)')],[0.5,0.5,0.5]);
                set(Fig.Cell_CDFLT_Err_a,'EdgeColor','None'); set(Fig.Cell_CDFLT_Err_a,'FaceAlpha',0.5);
                hold on; plot(CellValue.CDFBin'/(60/Param.Time.Exp),CellValue.CDFLT(:,1)','-k','LineWidth',1); hold off
                title('CDF: track length'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(CellValue.CDFLT(:,1)+CellValue.CDFLT(:,2)))        
            subplot(2,2,3)
                data = [CellValue.BI(1),CellValue.BM(1);CellValue.BI(3),CellValue.BM(3)]/(60/Param.Time.Exp);
                err = [CellValue.BI(2),CellValue.BM(2);CellValue.BI(4),CellValue.BM(4)]/(60/Param.Time.Exp);
                ngroups = size(data,1); nbars = size(data,2); groupwidth = min(0.8,nbars/(nbars+1.5));
                bar(data)
                hold on
                for i = 1:nbars
                    X = (1:ngroups)-groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
                    errorbar(X,data(:,i),err(:,i),'Color',[0,0,0],'LineStyle','None')
                end
                hold off
                set(gca,'xticklabel',{'Mean','Median'});
                legend({['IB, ',num2str(CellValue.BI(1)/(60/Param.Time.Exp),'%.2f'),'±',num2str(CellValue.BI(2)/(60/Param.Time.Exp),'%.2f'),' min'];...
                    ['MB, ',num2str(CellValue.BM(1)/(60/Param.Time.Exp),'%.2f'),'±',num2str(CellValue.BM(2)/(60/Param.Time.Exp),'%.2f'),' min']})
                title('Intrinsic and measured bursts'); ylabel('Burst width (min)'); ylim([-0.02,1.2]*max(max(data+err)))
            subplot(2,2,4)
                data = [CellValue.DI(1),CellValue.DM(1);CellValue.DI(3),CellValue.DM(3)]/(60/Param.Time.Exp);
                err = [CellValue.DI(2),CellValue.DM(2);CellValue.DI(4),CellValue.DM(4)]/(60/Param.Time.Exp);
                ngroups = size(data,1); nbars = size(data,2); groupwidth = min(0.8,nbars/(nbars+1.5));
                bar(data)
                hold on
                for i = 1:nbars
                    X = (1:ngroups)-groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
                    errorbar(X,data(:,i),err(:,i),'Color',[0,0,0],'LineStyle','None')
                end
                hold off
                set(gca,'xticklabel',{'Mean','Median'});
                legend({['ID, ',num2str(CellValue.DI(1)/(60/Param.Time.Exp),'%.2f'),'±',num2str(CellValue.DI(2)/(60/Param.Time.Exp),'%.2f'),' min'];...
                    ['MD, ',num2str(CellValue.DM(1)/(60/Param.Time.Exp),'%.2f'),'±',num2str(CellValue.DM(2)/(60/Param.Time.Exp),'%.2f'),' min']})
                title('Intrinsic and measured dwells'); ylabel('Durst width (min)'); if ~isnan(max(max(data+err))); ylim([-0.02,1.2]*max(max(data+err))); end
        Fig.tab9 = uitab(Fig.group, 'Title', 'SingleCell_2');
        Fig.axes9 = axes('Parent', Fig.tab9);
            subplot(2,2,1)
                Fig.Cell_CDFBM_Err_x = [CellValue.CDFBin'/(60/Param.Time.Exp),fliplr(CellValue.CDFBin'/(60/Param.Time.Exp))];
                Fig.Cell_CDFBM_Err_a = fill(Fig.Cell_CDFBM_Err_x,[CellValue.CDFBM(:,1)'+CellValue.CDFBM(:,2)',fliplr(CellValue.CDFBM(:,1)'-CellValue.CDFBM(:,2)')],[0.5,0.5,0.5]);
                set(Fig.Cell_CDFBM_Err_a,'EdgeColor','None'); set(Fig.Cell_CDFBM_Err_a,'FaceAlpha',0.5);
                hold on; plot(CellValue.CDFBin'/(60/Param.Time.Exp),CellValue.CDFBM(:,1)','-k','LineWidth',1); hold off
                legend(['FitMean: ',num2str(CellValue.BM(5)/(60/Param.Time.Exp),'%.1f'),'±',num2str(CellValue.BM(6)/(60/Param.Time.Exp),'%.1f'),' min'])
                title('CDF: measured burst'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(CellValue.CDFBM(:,1)+CellValue.CDFBM(:,2)))
            subplot(2,2,2)
                Fig.Cell_CDFDM_Err_x = [CellValue.CDFBin'/(60/Param.Time.Exp),fliplr(CellValue.CDFBin'/(60/Param.Time.Exp))];
                Fig.Cell_CDFDM_Err_a = fill(Fig.Cell_CDFDM_Err_x,[CellValue.CDFDM(:,1)'+CellValue.CDFDM(:,2)',fliplr(CellValue.CDFDM(:,1)'-CellValue.CDFDM(:,2)')],[0.5,0.5,0.5]);
                set(Fig.Cell_CDFDM_Err_a,'EdgeColor','None'); set(Fig.Cell_CDFDM_Err_a,'FaceAlpha',0.5);
                hold on; plot(CellValue.CDFBin'/(60/Param.Time.Exp),CellValue.CDFDM(:,1)','-k','LineWidth',1); hold off
                legend(['FitMean: ',num2str(CellValue.DM(5)/(60/Param.Time.Exp),'%.1f'),'±',num2str(CellValue.DM(6)/(60/Param.Time.Exp),'%.1f'),' min'])
                title('CDF: measured dwell'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)/6); ylim([-0.02,1.2]*max(CellValue.CDFDM(:,1)+CellValue.CDFDM(:,2)))
            subplot(2,2,3)
                Fig.Cell_CDFBMCor_Err_x = [CellValue.CDFBin'/(60/Param.Time.Exp),fliplr(CellValue.CDFBin'/(60/Param.Time.Exp))];
                Fig.Cell_CDFBMCor_Err_a = fill(Fig.Cell_CDFBMCor_Err_x,[CellValue.CDFBMCor(:,1)'+CellValue.CDFBMCor(:,2)',fliplr(CellValue.CDFBMCor(:,1)'-CellValue.CDFBMCor(:,2)')],[0.5,0.5,0.5]);
                set(Fig.Cell_CDFBMCor_Err_a,'EdgeColor','None'); set(Fig.Cell_CDFBMCor_Err_a,'FaceAlpha',0.5);
                hold on; plot(CellValue.CDFBin'/(60/Param.Time.Exp),CellValue.CDFBMCor(:,1)','-k','LineWidth',1); hold off
                legend(['CorFitMean: ',num2str(CellValue.BM(7)/(60/Param.Time.Exp),'%.1f'),'±',num2str(CellValue.BM(8)/(60/Param.Time.Exp),'%.1f'),' min'])
                title('CDF: apparent burst'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)); ylim([-0.02,1.2]*max(CellValue.CDFBMCor(:,1)+CellValue.CDFBMCor(:,2)))
            subplot(2,2,4)
                Fig.Cell_CDFDMCor_Err_x = [CellValue.CDFBin'/(60/Param.Time.Exp),fliplr(CellValue.CDFBin'/(60/Param.Time.Exp))];
                Fig.Cell_CDFDMCor_Err_a = fill(Fig.Cell_CDFDMCor_Err_x,[CellValue.CDFDMCor(:,1)'+CellValue.CDFDMCor(:,2)',fliplr(CellValue.CDFDMCor(:,1)'-CellValue.CDFDMCor(:,2)')],[0.5,0.5,0.5]);
                set(Fig.Cell_CDFDMCor_Err_a,'EdgeColor','None'); set(Fig.Cell_CDFDMCor_Err_a,'FaceAlpha',0.5);
                hold on; plot(CellValue.CDFBin'/(60/Param.Time.Exp),CellValue.CDFDMCor(:,1)','-k','LineWidth',1); hold off
                legend(['CorFitMean: ',num2str(CellValue.DM(7)/(60/Param.Time.Exp),'%.1f'),'±',num2str(CellValue.DM(8)/(60/Param.Time.Exp),'%.1f'),' min'])
                title('CDF: apparent dwell'); xlabel('Time (min)'); ylabel('Frequency'); xlim([-0.02,1.02]*Param.Frame.LTmax/(60/Param.Time.Exp)/6); ylim([-0.02,1.2]*max(CellValue.CDFDMCor(:,1)+CellValue.CDFDMCor(:,2)))
    end
    clear ans data err dmy1 dmy2 ngroups nbars groupwidth X
    if Mode.ISP == 1
        for i = 1:N
            disp(['i = ',num2str(i)])
            Fig.axes1 = axes('Parent', Fig.tab1);
            subplot(3,1,1)
                plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).nR,'Color','#666666','LineWidth',1)
                hold on; plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).TSI*max(Data(i).nR),'-k','LineWidth',1); hold off
                title('ITS (black) and #R (gray)'); ylabel('#R'); xlim([0,ceil(Param.Frame.LTmax/(60/Param.Time.Exp))]); ylim([-0.02,1.2]*max(Data(i).nR))
            subplot(3,1,2)
                plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).nT1,'Color','#666666','LineWidth',1)
                hold on; plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).nT2,'Color','#999999','LineWidth',1); hold off
                title('#T1 (gray) and #T2 (light gray)'); ylabel('#Tags'); xlim([0,ceil(Param.Frame.LTmax/(60/Param.Time.Exp))]); ylim([-0.02,1.2]*max(Data(i).nT1))
            subplot(3,1,3)
                plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).FIT1,'Color','#666666','LineWidth',1)
                hold on
                plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).FIT2,'Color','#999999','LineWidth',1)
                plot((1:Data(i).LT)/(60/Param.Time.Exp),Data(i).FIT1B*max(Data(i).FIT1),'Color','#000000','LineWidth',1)
                if Mode.FRAP == 1
                    plot([1,Data(i).LT-Param.Frame.FRAP-1]/(60/Param.Time.Exp),[1,1]*Value.FI_Mean_Theory,':r','LineWidth',1);
                    plot((Data(i).LT-Param.Frame.FRAP:1:Data(i).LT)/(60/Param.Time.Exp),FRAP.Theory*Value.FI_Mean_Theory,'--r','LineWidth',1);
                end
                if Mode.RO == 1
                    plot([1,Data(i).LT-Param.Frame.RO-1]/(60/Param.Time.Exp),[1,1]*Value.FI_Mean_Theory,':r','LineWidth',1);
                    plot((Data(i).LT-Param.Frame.RO:1:Data(i).LT)/(60/Param.Time.Exp),RO.Theory*Value.FI_Mean_Theory,'--r','LineWidth',1);
                end
                hold off
                title('FIT1 (gray), FIT2 (light gray), and MTS (black)'); xlabel('Time (min)'); ylabel('FI (AU)'); xlim([0,ceil(Param.Frame.LTmax/(60/Param.Time.Exp))]); ylim([-0.02,1.2]*max(Data(i).FIT1))
                pause
        end
    end
end
clear i

%% Functions
function Data = runmc(Mode,Param,L)
    % Decompose variables for parallel computation
    modeE = Mode.E; modeB = Mode.B; modeLT = Mode.LT; modeFRAP = Mode.FRAP; modeRO = Mode.RO; modenR0 = Mode.nR0; modenRC = Mode.nRC;
    ki = Param.Rate.i; ke = Param.Rate.e; keA = Param.Rate.eA; kt = Param.Rate.t; kon = Param.Rate.on; koff = Param.Rate.off;
    LTmin = Param.Frame.LTmin; LTmu = Param.Frame.LTmu; LTmax = Param.Frame.LTmax; fEq = Param.Frame.Eq; fFRAP = Param.Frame.FRAP; fRO = Param.Frame.RO; fCO = Param.Frame.CO;
    tdT = Param.Time.dT; tExp = Param.Time.Exp;
    l1 = Param.RNA.l1; nt1 = Param.RNA.nt1; lt1 = Param.RNA.lt1; l2 = Param.RNA.l2; lpA = Param.RNA.lpA; nt2 = Param.RNA.nt2; lt2 = Param.RNA.lt2; l3 = Param.RNA.l3;
    pShape = Param.Photon.shape; pScale = Param.Photon.scale; pTag = Param.Photon.tag;
    RFP = Param.RFP; BTS = Param.BTS; ITS = Param.ITS;
    %--- Codon-based effective ke and #Tags
    c_ke = ke*ones(L,1); % Elongation rate
    c_ke(l1+nt1*lt1+l2+1:l1+nt1*lt1+l2+lpA) = keA;
    c_nT1 = (1:1:L)'-l1-1; % Tag1
    c_nT1(1:l1) = 0;
    c_nT1(l1+1:l1+nt1*lt1) = floor(c_nT1(l1+1:l1+nt1*lt1)/lt1);
    c_nT1(l1+nt1*lt1+1:end) = nt1;
    n = l1+nt1*lt1+l2+lpA;
    c_nT2 = (1:1:L)'-n-1; % Tag2
    c_nT2(1:n) = 0;
    c_nT2(n+1:n+nt2*lt2) = floor(c_nT2(n+1:n+nt2*lt2)/lt2);
    c_nT2(n+nt2*lt2+1:end) = nt2;
    %--- Determine LT and simulation steps
    if modeLT == 0; LT = LTmax;
    else; LT = min(LTmax,round(exprnd(LTmu)+LTmin));
    end
    nS = (fEq+LT)*tExp/tdT;
    %--- Initialize variables
    TSI = nan(1,(fEq+LT)*tExp/tdT); % Intrinsic translational status
    nR = nan(1,(fEq+LT)*tExp/tdT); % # Ribosome
    nT1 = nan(1,(fEq+LT)*tExp/tdT); % #Tag1
    nT2 = nan(1,(fEq+LT)*tExp/tdT); % #Tag2
    FIT1 = nan(1,(fEq+LT)*tExp/tdT); % Fluorescence intensity of Tag1
    FIT2 = nan(1,(fEq+LT)*tExp/tdT); % Fluorescence intensity of Tag2
    CD = 0; % Collision density
    CO = zeros(1,L); % Codon occupancy
    iR = []; % Ribosome information of each RNA
             % [Index,Initiated frame,Location,...
             % #Active Tag1,#Total pre-Tag1,#Total Tag1,...
             % #Active Tag2,#Total pre-Tag2,#Total Tag2,...
    iMP = []; % Mature protein information of each RNA
              % [iR;Termination case;Terminated frame]
    idxNP = 0; % Index of nascent peptide
    BCD = nan(1,(fEq+LT)*tExp/tdT); % Bursting countdown
    ICD = nan(1,(fEq+LT)*tExp/tdT); % Initiation countdown
    udt = []; % Ribosome update order
    %--- Simulation at each simulation step; 'j' = simulation step
    for j = 1:nS
        %--- Initialize parameters
        nT1(j) = 0; nT2(j) = 0;
        %--- Countdowns
        if j>1
            BCD(j) = BCD(j-1)-tdT;
            ICD(j) = ICD(j-1)-tdT;
        end
        %--- Determine intrinsic translational status
        if modeB == 0 % Constitutive
            TSI(j) = 1;
            if j == 1; ICD(j) = gamrnd(ITS,1/ki/ITS); end
        elseif modeB == 1 % 2SRTwOV
            if j == 1
                if rand <= kon/(kon+koff); TSI(j) = 1; BCD(j) = gamrnd(BTS,1/koff/BTS); ICD(j) = gamrnd(ITS,1/ki/ITS);
                else; TSI(j) = 0; 
                end
            else
                if TSI(j-1) == 1 && BCD(j) <= 0; TSI(j) = 0; BCD(j) = nan; ICD(j) = nan;
                elseif TSI(j-1) == 0 && rand <= 1-exp(-kon*tdT); TSI(j) = 1; BCD(j) = gamrnd(BTS,1/koff/BTS); ICD(j) = gamrnd(ITS,1/ki/ITS);
                else; TSI(j) = TSI(j-1);
                end
            end
        elseif modeB == 2 % 2SRTwoOV
            if j == 1
                if rand <= kon/(kon+koff); TSI(j) = 1; BCD(j) = gamrnd(BTS,1/koff/BTS); ICD(j) = gamrnd(ITS,1/ki/ITS);
                else; TSI(j) = 0; 
                end
            else
                if TSI(j-1) == 1 && BCD(j) <= 0; TSI(j) = 0; BCD(j) = nan; ICD(j) = nan;
                elseif TSI(j-1) == 0 && nR(j-1) == 0 && rand <= 1-exp(-kon*tdT); TSI(j) = 1; BCD(j) = gamrnd(BTS,1/koff/BTS); ICD(j) = gamrnd(ITS,1/ki/ITS);
                else; TSI(j) = TSI(j-1);
                end
            end
        end            
        if modeRO == 1 && j >= (fEq+LT-(fRO+1))*tExp/tdT+1 % Off-state for RO
            TSI(j) = 0; ICD(j) = nan;
        end
        %--- Initiation
        if j == 1; nR(j) = 0; % No initiation at first step
        else
            if ICD(j) <= 0
                if nR(j-1) == 0 || (~isempty(iR) && iR(end,3) > RFP)
                    nR(j) = nR(j-1)+1; idxNP = idxNP+1;
                    iR(nR(j),1) = idxNP; iR(nR(j),2) = j; iR(nR(j),3) = 1;
                    ICD(j) = gamrnd(ITS,1/ki/ITS);
                elseif ~isempty(iR) && iR(end,3) <= RFP
                    nR(j) = nR(j-1); CD = CD+1;
                end
            else; nR(j) = nR(j-1);
            end
        end
        %--- Elongation
        if nR(j) ~= 0
            %--- Determine update order
            if modeE == 0 || modeE == 2; udt = (1:nR(j));
            elseif modeE == 1 || modeE == 3; udt = randperm(nR(j),nR(j));
            end
            %-- TC movement and #Tags
            for k = 1:nR(j)
                idx = find(udt==k);
                %--- TC movement
                if modeE == 0 || modeE == 1 % Deterministic elongation
                    if iR(idx,2) ~= j % Ribosomes that are not initiated at the same step
                        if idx == 1 && iR(idx,3) < L % Leading ribosome that is not at the end
                            iR(idx,3) = iR(idx,3)+c_ke(round(iR(idx,3)))*tdT;
                        elseif idx > 1 && iR(idx-1,3)-(iR(idx,3)+c_ke(round(iR(idx,3)))*tdT) < RFP % Collision
                            CD = CD+1;
                        elseif idx > 1 && iR(idx-1,3)-(iR(idx,3)+c_ke(round(iR(idx,3)))*tdT) >= RFP % Movement
                            iR(idx,3) = iR(idx,3)+c_ke(round(iR(idx,3)))*tdT;
                        end
                    end
                elseif modeE == 2 || modeE == 3 % TASEP
                    if iR(idx,2) ~= j
                        if idx == 1 && iR(idx,3) < L && rand <= 1-exp(-c_ke(round(iR(idx,3)))*tdT)
                            iR(idx,3) = iR(idx,3)+1;
                        elseif idx > 1 && iR(idx-1,3)-(iR(idx,3)+1) < RFP
                            CD = CD+1;
                        elseif idx > 1 && iR(idx-1,3)-iR(idx,3) >= RFP && rand <= 1-exp(-c_ke(round(iR(idx,3)))*tdT)
                            iR(idx,3) = iR(idx,3)+1;                           
                        end
                    end
                end
                if iR(idx,3) > L; iR(idx,3) = L; end
                %--- #Tags
                if round(iR(idx,3)) <= l1+lt1; iR(idx,4) = 0; iR(idx,5) = 0; iR(idx,6) = 0;
                else; iR(idx,5) = iR(idx,6); iR(idx,6) = c_nT1(round(iR(idx,3)));
                end
                iR(idx,4) = iR(idx,4)+iR(idx,6)-iR(idx,5);
                nT1(j) = nT1(j)+iR(idx,4);
                if round(iR(idx,3)) <= (l1+nt1*lt1+l2+lpA+nt2*lt2); iR(idx,7) = 0; iR(idx,8) = 0; iR(idx,9) = 0;
                else; iR(idx,8) = iR(idx,9); iR(idx,9) = c_nT2(round(iR(idx,3)));
                end
                iR(idx,7) = iR(idx,7)+iR(idx,9)-iR(idx,8);
                nT2(j) = nT2(j)+iR(idx,7);
                %--- Codon occupancy
                CO(round(iR(idx,3))) = CO(round(iR(idx,3)))+1;
            end
        end
        %--- Termination
        if nR(j) ~= 0 && iR(1,3) == L && rand <= 1-exp(-kt*tdT)
            nR(j) = nR(j)-1; nT1(j) = nT1(j)-iR(1,4); nT2(j) = nT2(j)-iR(1,7);
            iMP = [iMP;[iR(1,:),0,j]];
            iR(1,:) = [];
        end
        %--- Fluorescence intensity
        FIT1(j) = gamrnd(nT1(j)*pShape,pScale);
        FIT2(j) = gamrnd(nT2(j)*pShape,pScale);
        %--- FRAP
        if modeFRAP == 1 && j == (fEq+LT-(fFRAP+1))*tExp/tdT+1
            for k = 1:nR(j); iR(k,4) = 0; iR(k,7) = 0; end
            FIT1(j) = 0; FIT2(j) = 0;
        end
    end
    %--- Reshape time-dependent informations
    TSI = reshape(TSI,[tExp/tdT,fEq+LT])'; TSI = TSI(fEq+1:fEq+LT,1);
    ICD = reshape(ICD,[tExp/tdT,fEq+LT])'; ICD = ICD(fEq+1:fEq+LT,1);
    nR = reshape(nR,[tExp/tdT,fEq+LT])'; nR = nR(fEq+1:fEq+LT,1);
    nT1 = reshape(nT1,[tExp/tdT,fEq+LT])'; nT1 = nT1(fEq+1:fEq+LT,1);
    nT2 = reshape(nT2,[tExp/tdT,fEq+LT])'; nT2 = nT2(fEq+1:fEq+LT,1);
    FIT1 = reshape(FIT1,[tExp/tdT,fEq+LT])'; FIT1 = FIT1(fEq+1:fEq+LT,1);
    FIT2 = reshape(FIT2,[tExp/tdT,fEq+LT])'; FIT2 = FIT2(fEq+1:fEq+LT,1);
    %--- #R and FIT1 sampling for every 'gapSample'-frames
    gapSample = 60; % In 'frame'
    fSample = (1:gapSample:floor((LT-1)/gapSample)*gapSample+1);
    if modeFRAP == 1; fSample = (1:gapSample:floor(((LT-fFRAP)-1)/gapSample)*gapSample+1); end
    if modeRO == 1; fSample = (1:gapSample:floor(((LT-fRO)-1)/gapSample)*gapSample+1); end
    nRS = nR(fSample); FIT1S = FIT1(fSample);
    %--- Normalize CD and CO
    CD = tExp*CD/nS/L/tdT; % #Collisions per minute per frame
    CO = CO/sum(CO); % Normalized to 1
    %--- Filtering and binarization of FIT1
    FIT1F = medfilt1(FIT1);
    FIT1B = zeros(LT,1); FIT1B(FIT1F>=0.05*max(FIT1F)) = 1;
    %--- Ignore shorter measured bursts than fCO frames
    idx = find(FIT1B(1:end-1)~=FIT1B(2:end));
    if ~isempty(idx)
        if idx(1) <= fCO; idx = [0;idx]; elseif LT-idx(end) <= fCO; idx = [idx;LT]; end
        idx2 = find(idx(2:end)-idx(1:end-1)<=fCO)+1;
        if ~isempty(idx2)
            for k = 1:length(idx2)
                if FIT1B(idx(idx2(k))) == 1
                    FIT1B(idx(idx2(k)-1)+1:idx(idx2(k))) = ~FIT1B(idx(idx2(k)-1)+1:idx(idx2(k)));
                end
            end
        end
    end
    %--- Ignore shorter measured dwells than fCO frames
    idx = find(FIT1B(1:end-1)~=FIT1B(2:end));
    if ~isempty(idx)
        if idx(1) <= fCO; idx = [0;idx]; elseif LT-idx(end) <= fCO; idx = [idx;LT]; end
        idx2 = find(idx(2:end)-idx(1:end-1)<=fCO)+1;
        if ~isempty(idx2)
            for k = 1:length(idx2)
                if FIT1B(idx(idx2(k))) == 0
                    FIT1B(idx(idx2(k)-1)+1:idx(idx2(k))) = ~FIT1B(idx(idx2(k)-1)+1:idx(idx2(k)));
                end
            end
        end
    end
    %--- Get burst characteristics
    dI = diff(TSI); % Intrinsic
    idx = find(dI~=0);
    if isempty(idx)
        BI = [1,LT,TSI(1),LT]; % [Start,End,Status,Length]
    else
        BI = [[1;idx+1],[idx;LT],[TSI(1);dI(idx)]]; BI(:,4) = BI(:,2)-BI(:,1)+1;
        BI(BI(:,3)==-1,3) = 0;
    end
    dM = diff(FIT1B); % Measured
    idx = find(dM~=0);
    if isempty(idx)
        BM = [1,LT,FIT1B(1),LT,mean(FIT1),max(FIT1)]; % [Start,End,Status,Length,MeanFI,MaxFI]
    else
        BM = [[1;idx+1],[idx;LT],[FIT1B(1);dM(idx)]]; BM(:,4) = BM(:,2)-BM(:,1)+1;
        for k = 1:length(idx)+1
            BM(k,5) = mean(FIT1(BM(k,1):BM(k,2))); BM(k,6) = max(FIT1(BM(k,1):BM(k,2)));
        end
        BM(BM(:,3)==-1,3) = 0;
    end
    %--- Summarize simulation result into a structure
    Data.LT = LT; % Track length
    Data.nS = nS; % Simulation step
    Data.TSI = TSI; % Intrinsic translational status
    Data.ICD = ICD; % History of initiation events
    Data.nR = nR; % #Ribosome
    Data.nT1 = nT1; % #Tag1
    Data.nT2 = nT2; % #Tag2
    Data.FIT1 = FIT1; % FI Tag1
    Data.FIT1B = FIT1B; % Measured translational status
    Data.FIT2 = FIT2; % FI Tag2
    Data.CD = CD; % Collision density
    Data.CO = CO; % Codon occupancy
    Data.iR = iR; % Ribosome information
    Data.iMP = iMP; % Mature protein information
    Data.nRS = nRS; % Sampled #Ribosome
    Data.FIT1S = FIT1S; % Sampled FIT1 for #R
    Data.BI = BI; % Intrinsic burst characteristics
    Data.BM = BM; % Measured burst characteristics
end