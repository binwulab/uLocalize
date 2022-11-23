function gfit = pdFitGlobal(arginCI,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp)

    vn_fit = {'MBExp','MDExp','NPExp','MBFit','MDFit','NPFit','ki','koff','kon','init/min','IB','ID','idxki','idxIB','idxID','RMSE'};
    vn_ci = {'ki_LB','ki_UB','koff_LB','koff_UB','kon_LB','kon_UB','init/min_LB','init/min_UB','IB_LB','IB_UB','ID_LB','ID_UB',};

    % ORF series
    idxexp = (1:4);
    fit = zeros(length(idxexp),16);
    nrmse = sum(PDnRMSE(:,:,:,idxexp,4),4);
    [var,idx] = min(nrmse,[],'all','linear');   
    [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);

    for i = 1:length(idxexp)
        fit(i,:) = [arrayExp(idxexp(i),[1,3,5]),squeeze(PD(i1,i2,i3,arrayL(idxexp(i)),[1,3,5]))',IN.ki(i1),IN.koff(i2),IN.kon(i3),...
            IN.nInit(i1),IN.IB(i2),IN.ID(i3),i1,i2,i3,PDnRMSE(i1,i2,i3,idxexp(i),4)];
    end

    gfit.ORF = array2table(fit,"VariableNames",vn_fit);
    gfit.ORF.Properties.RowNames = resExp.Properties.RowNames(idxexp);
    
    ci = zeros(length(idxexp),12);

    for i = 1:length(idxexp)
        i1 = find(abs(nrmse(:,fit(i,14),fit(i,15))-var)<=arginCI*length(idxexp));
        i2 = find(abs(nrmse(fit(i,13),:,fit(i,15))-var)<=arginCI*length(idxexp));
        i3 = find(abs(nrmse(fit(i,13),fit(i,14),:)-var)<=arginCI*length(idxexp));
        ci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
            IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
    end

    gfit.ORF_CI = array2table(ci,"VariableNames",vn_ci);
    gfit.ORF_CI.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    % Hairpin series, Model 1 (HP affects the burst timing)
    idxexp = [1,5:7];
    fit = zeros(length(idxexp),16);
    var = zeros(N.ki,length(idxexp));
    idx = zeros(N.ki,length(idxexp));

    for i = 1:N.ki
        for j = 1:length(idxexp)
            [var(i,j),idx(i,j)] = min(PDnRMSE(i,:,:,idxexp(j),4),[],'all','linear');
        end
    end

    varsum = sum(var,2);
    [var2,i1] = min(varsum);
    idx2 = idx(i1,:)';
    [i2,i3] = ind2sub([N.koff,N.kon],idx2);

    for i = 1:length(idxexp)
        fit(i,:) = [arrayExp(idxexp(i),[1,3,5]),squeeze(PD(i1,i2(i),i3(i),arrayL(idxexp(i)),[1,3,5]))',IN.ki(i1),IN.koff(i2(i)),IN.kon(i3(i)),...
            IN.nInit(i1),IN.IB(i2(i)),IN.ID(i3(i)),i1,i2(i),i3(i),PDnRMSE(i1,i2(i),i3(i),idxexp(i),4)];
    end

    gfit.HP_M1 = array2table(fit,"VariableNames",vn_fit);
    gfit.HP_M1.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    ci = zeros(length(idxexp),12);

    for i = 1:length(idxexp)
        i1 = find(abs(varsum-var2)<=arginCI*length(idxexp));
        i2 = find(abs(PDnRMSE(fit(i,13),:,fit(i,15),idxexp(i),4)-var(fit(i,13),i))<=arginCI);
        i3 = find(abs(PDnRMSE(fit(i,13),fit(i,14),:,idxexp(i),4)-var(fit(i,13),i))<=arginCI);
        ci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
            IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
    end

    gfit.HP_M1_CI = array2table(ci,"VariableNames",vn_ci);
    gfit.HP_M1_CI.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    % Hairpin series, Model 2 (HP affects the burst amplitude)
    idxexp = [1,5:7];
    fit = zeros(length(idxexp),16);
    var = zeros(N.koff,N.kon,length(idxexp));
    idx = zeros(N.koff,N.kon,length(idxexp));

    for i = 1:N.koff
        for j = 1:N.kon
            for k = 1:length(idxexp)
                [var(i,j,k),idx(i,j,k)] = min(PDnRMSE(:,i,j,idxexp(k),4));
            end
        end
    end

    varsum = sum(var,3);
    [var2,idx2] = min(varsum,[],'all','linear');
    [i2,i3] = ind2sub([N.koff,N.kon],idx2);
    i1 = squeeze(idx(i2,i3,:));

    for i = 1:length(idxexp)
        fit(i,:) = [arrayExp(idxexp(i),[1,3,5]),squeeze(PD(i1(i),i2,i3,arrayL(idxexp(i)),[1,3,5]))',IN.ki(i1(i)),IN.koff(i2),IN.kon(i3),...
            IN.nInit(i1(i)),IN.IB(i2),IN.ID(i3),i1(i),i2,i3,PDnRMSE(i1(i),i2,i3,idxexp(i),4)];
    end

    gfit.HP_M2 = array2table(fit,"VariableNames",vn_fit);
    gfit.HP_M2.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    ci = zeros(length(idxexp),12);

    for i = 1:length(idxexp)
        i1 = find(abs(PDnRMSE(:,fit(i,14),fit(i,15),idxexp(i),4)-var(fit(i,14),fit(i,15),i))<=arginCI);
        i2 = find(abs(varsum(:,fit(i,15))-var2)<=arginCI*length(idxexp));
        i3 = find(abs(varsum(fit(i,14),:)-var2)<=arginCI*length(idxexp));
        ci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
            IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
    end

    gfit.HP_M2_CI = array2table(ci,"VariableNames",vn_ci);
    gfit.HP_M2_CI.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    % Torin series, ST-AID
    arrayExp(1,[9,10]) = [9.460,1.064]; % Use FISH results with vesicle treatment
    PDnRMSE = pdGetRMSE(N,PD,arrayExp,arrayL); % Recalculate PDnRMSE

    idxexp = [1,13];
    fit = zeros(length(idxexp),16);
    var = zeros(N.koff,length(idxexp));
    idx = zeros(N.koff,length(idxexp));

    for i = 1:N.koff
        for j = 1:length(idxexp)
            [var(i,j),idx(i,j)] = min(PDnRMSE(:,i,:,idxexp(j),4),[],'all','linear');
        end
    end

    varsum = sum(var,2);
    [var2,i2] = min(varsum);
    idx2 = idx(i2,:)';
    [i1,i3] = ind2sub([N.ki,N.kon],idx2);

    for i = 1:length(idxexp)
        fit(i,:) = [arrayExp(idxexp(i),[1,3,5]),squeeze(PD(i1(i),i2,i3(i),arrayL(idxexp(i)),[1,3,5]))',IN.ki(i1(i)),IN.koff(i2),IN.kon(i3(i)),...
            IN.nInit(i1(i)),IN.IB(i2),IN.ID(i3(i)),i1(i),i2,i3(i),PDnRMSE(i1(i),i2,i3(i),idxexp(i),4)];
    end

    staidfit = array2table(fit,"VariableNames",vn_fit);
    staidfit.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    ci = zeros(length(idxexp),12);

    for i = 1:length(idxexp)
        i1 = find(abs(PDnRMSE(:,fit(i,14),fit(i,15),idxexp(i),4)-var(fit(i,14),i))<=arginCI);
        i2 = find(abs(varsum-var2)<=arginCI*length(idxexp));
        i3 = find(abs(PDnRMSE(fit(i,13),fit(i,14),:,idxexp(i),4)-var(fit(i,14),i))<=arginCI);
        ci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
            IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
    end

    staidci = array2table(ci,"VariableNames",vn_ci);
    staidci.Properties.RowNames = resExp.Properties.RowNames(idxexp);    

    % Torin series, EEF1A
    idxexp = [14,15];
    fit = zeros(length(idxexp),16);
    var = zeros(N.koff,length(idxexp));
    idx = zeros(N.koff,length(idxexp));

    for i = 1:N.koff
        for j = 1:length(idxexp)
            [var(i,j),idx(i,j)] = min(PDnRMSE(:,i,:,idxexp(j),4),[],'all','linear');
        end
    end

    varsum = sum(var,2);
    [var2,i2] = min(varsum);
    idx2 = idx(i2,:)';
    [i1,i3] = ind2sub([N.ki,N.kon],idx2);

    for i = 1:length(idxexp)
        fit(i,:) = [arrayExp(idxexp(i),[1,3,5]),squeeze(PD(i1(i),i2,i3(i),arrayL(idxexp(i)),[1,3,5]))',IN.ki(i1(i)),IN.koff(i2),IN.kon(i3(i)),...
            IN.nInit(i1(i)),IN.IB(i2),IN.ID(i3(i)),i1(i),i2,i3(i),PDnRMSE(i1(i),i2,i3(i),idxexp(i),4)];
    end

    eef1afit = array2table(fit,"VariableNames",vn_fit);
    eef1afit.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    ci = zeros(length(idxexp),12);

    for i = 1:length(idxexp)
        i1 = find(abs(PDnRMSE(:,fit(i,14),fit(i,15),idxexp(i),4)-var(fit(i,14),i))<=arginCI);
        i2 = find(abs(varsum-var2)<=arginCI*length(idxexp));
        i3 = find(abs(PDnRMSE(fit(i,13),fit(i,14),:,idxexp(i),4)-var(fit(i,14),i))<=arginCI);
        ci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
            IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
    end
    
    eef1aci = array2table(ci,"VariableNames",vn_ci);
    eef1aci.Properties.RowNames = resExp.Properties.RowNames(idxexp); 

    % Summarize results of Torin series
    gfit.Torin = [staidfit;eef1afit];
    gfit.Torin_CI = [staidci;eef1aci];

end

