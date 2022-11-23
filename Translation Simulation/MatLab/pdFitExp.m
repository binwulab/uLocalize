function resfit = pdFitExp(arginCI,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp)

    vn_fit = {'MBExp','MDExp','NPExp','MBFit','MDFit','NPFit','ki','koff','kon','init/min','IB','ID','idxki','idxIB','idxID','RMSE'};
    vn_ci = {'ki_LB','ki_UB','koff_LB','koff_UB','kon_LB','kon_UB','init/min_LB','init/min_UB','IB_LB','IB_UB','ID_LB','ID_UB'};

    fit = zeros(size(arrayExp,1),16);
    ci = zeros(size(arrayExp,1),12);

    for i = 1:size(arrayExp,1)
        
        % Get global minimum
        [var,idx] = min(PDnRMSE(:,:,:,i,4),[],'all','linear');   
        [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);

        fit(i,:) = [arrayExp(i,[1,3,5]),squeeze(PD(i1,i2,i3,arrayL(i),[1,3,5]))',IN.ki(i1),IN.koff(i2),IN.kon(i3),...
            IN.nInit(i1),IN.IB(i2),IN.ID(i3),i1,i2,i3,PDnRMSE(i1,i2,i3,i,4)];
        
        % Calculate CI
        i1 = find(abs(PDnRMSE(:,fit(i,14),fit(i,15),i,4)-var)<=arginCI);
        i2 = find(abs(PDnRMSE(fit(i,13),:,fit(i,15),i,4)-var)<=arginCI);
        i3 = find(abs(PDnRMSE(fit(i,13),fit(i,14),:,i,4)-var)<=arginCI);


        ci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
                IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];

    end

    resfit.Fit = array2table(fit,"VariableNames",vn_fit);
    resfit.Fit.Properties.RowNames = resExp.Properties.RowNames;

    resfit.CI = array2table(ci,"VariableNames",vn_ci);
    resfit.CI.Properties.RowNames = resExp.Properties.RowNames;

end