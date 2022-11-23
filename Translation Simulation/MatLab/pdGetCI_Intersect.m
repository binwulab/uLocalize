function [resfit,resgfit] = pdGetCI_Intersect(arginCI_Intersect,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp,resFit,resGFit)
    
    resfit = resFit;
    resgfit = resGFit;

    vn_idxci = {'idxMB','idxMD','idxNP','idxIntersect'};
    vn_resci = {'ki_LB','ki_UB','koff_LB','koff_UB','kon_LB','kon_UB','init/min_LB','init/min_UB','IB_LB','IB_UB','ID_LB','ID_UB',};
    
    if strcmp(arginCI_Intersect{1},'SEM')
        nrmse = PDnRMSE;
    elseif strcmp(arginCI_Intersect{1},'Mean') % Normalize PDRMSE with the mean values
        nrmse = zeros(N.ki,N.koff,N.kon,size(arrayExp,1),4);
        for i = 1:size(arrayExp,1)
            for j = 1:3
                wgt = arrayExp(i,2*(j-1)+1);
                nrmse(:,:,:,i,j) = (sqrt((PD(:,:,:,arrayL(i),2*(j-1)+1)-arrayExp(i,2*(j-1)+1)).^2))./wgt;
            end
            nrmse(:,:,:,i,4) = nrmse(:,:,:,i,1)+nrmse(:,:,:,i,2)+nrmse(:,:,:,i,3);
        end
    end

    thr = arginCI_Intersect{2};

    % CI of individual fit results
    idxci = cell(size(arrayExp,1),4);
    resci = zeros(size(arrayExp,1),12);
    for i = 1:size(arrayExp,1)
        for j = 1:3 % CI_MB, CI_MD, and CI_NP
            idx = find(nrmse(:,:,:,i,j)<=thr);
            [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);
            idxci{i,j} = [idx,i1,i2,i3];
        end
        % CI_Overlap
        idx = intersect(intersect(idxci{i,1}(:,1),idxci{i,2}(:,1)),idxci{i,3}(:,1));
        [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);
        idxci{i,4} = [idx,i1,i2,i3];
        if isempty(idx)
            resci(i,:) = nan(1,12);
        else
            resci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
                IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
        end
    end
    resfit.CI_Intersect.Index = cell2table(idxci,"VariableNames",vn_idxci);
    resfit.CI_Intersect.Index.Properties.RowNames = resExp.Properties.RowNames;
    resfit.CI_Intersect.CI = array2table(resci,"VariableNames",vn_resci);
    resfit.CI_Intersect.CI.Properties.RowNames = resExp.Properties.RowNames;
    
    % CI of global fit results - ORF series
    idxexp = (1:4);
    idx = intersect(intersect(idxci{idxexp(1),4}(:,1),idxci{idxexp(2),4}(:,1)),intersect(idxci{idxexp(3),4}(:,1),idxci{idxexp(4),4}(:,1)));
    [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);
    if isempty(idx)
        resci = nan(1,12);
    else
        resci = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
                    IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
    end
    resgfit.CI_Intersect.ORF_Index = cell2table({[idx,i1,i2,i3]},"VariableNames",{'idxOverlap'});
    resgfit.CI_Intersect.ORF_CI = array2table(resci,"VariableNames",vn_resci);


    % CI of global fit results - HP series, Model 1
    idxexp = [1,5:7];
    idx_shared = intersect(intersect(idxci{idxexp(1),3}(:,1),idxci{idxexp(2),3}(:,1)),intersect(idxci{idxexp(3),3}(:,1),idxci{idxexp(4),3}(:,1)));
    idxhp = cell(4,1);
    resci = zeros(size(idxexp,1),12);
    for i = 1:length(idxexp)
        idx = intersect(intersect(idxci{idxexp(i),1},idxci{idxexp(i),2}),idx_shared);
        [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);
        idxhp{i,1} = [idx,i1,i2,i3];
        if isempty(idx)
            resci(i,:) = nan(1,12);
        else
            resci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
                IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
        end
    end
    resgfit.CI_Intersect.HP_M1_Index = cell2table(idxhp,"VariableNames",{'idxOverlap'});
    resgfit.CI_Intersect.HP_M1_Index.Properties.RowNames = resExp.Properties.RowNames(idxexp);
    resgfit.CI_Intersect.HP_M1_CI = array2table(resci,"VariableNames",vn_resci);
    resgfit.CI_Intersect.HP_M1_CI.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    % CI of global fit results - HP series, Model 2
    idx_shared = intersect(intersect(intersect(idxci{idxexp(1),1}(:,1),idxci{idxexp(2),1}(:,1)),intersect(idxci{idxexp(3),1}(:,1),idxci{idxexp(4),1}(:,1))),...
        intersect(intersect(idxci{idxexp(1),2}(:,1),idxci{idxexp(2),2}(:,1)),intersect(idxci{idxexp(3),2}(:,1),idxci{idxexp(4),2}(:,1))));
    for i = 1:length(idxexp)
        idx = intersect(idx_shared,idxci{idxexp(i),3}(:,1));
        [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);
        idxhp{i,1} = [idx,i1,i2,i3];
        if isempty(idx)
            resci(i,:) = nan(1,12);
        else
            resci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
                IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
        end
    end
    resgfit.CI_Intersect.HP_M2_Index = cell2table(idxhp,"VariableNames",{'idxOverlap'});
    resgfit.CI_Intersect.HP_M2_Index.Properties.RowNames = resExp.Properties.RowNames(idxexp);
    resgfit.CI_Intersect.HP_M2_CI = array2table(resci,"VariableNames",vn_resci);
    resgfit.CI_Intersect.HP_M2_CI.Properties.RowNames = resExp.Properties.RowNames(idxexp);

    % CI of global fit results - Torin series
    idxexp = [1,13:15];
    idxtorin = cell(4,1);
    idx_shared_staid = intersect(idxci{idxexp(1),2}(:,1),idxci{idxexp(2),3}(:,1));
    idx_shared_eef1a = intersect(idxci{idxexp(3),2}(:,1),idxci{idxexp(4),3}(:,1));
    for i = 1:length(idxexp)
        if i <= 2
            idx_shared = idx_shared_staid;
        elseif i >= 3
            idx_shared = idx_shared_eef1a;
        end
        idx = intersect(idx_shared,intersect(idxci{idxexp(i),1}(:,1),idxci{idxexp(i),3}(:,1)));
        [i1,i2,i3] = ind2sub([N.ki,N.koff,N.kon],idx);
        idxtorin{i,1} = [idx,i1,i2,i3];
        if isempty(idx)
            resci(i,:) = nan(1,12);
        else
            resci(i,:) = [IN.ki([min(i1),max(i1)])',IN.koff([min(i2),max(i2)])',IN.kon([min(i3),max(i3)])',...
                IN.nInit([min(i1),max(i1)])',IN.IB([max(i2),min(i2)])',IN.ID([max(i3),min(i3)])'];
        end
    end
    resgfit.CI_Intersect.Torin_Index = cell2table(idxtorin,"VariableNames",{'idxOverlap'});
    resgfit.CI_Intersect.Torin_Index.Properties.RowNames = resExp.Properties.RowNames(idxexp);
    resgfit.CI_Intersect.Torin_CI = array2table(resci,"VariableNames",vn_resci);
    resgfit.CI_Intersect.Torin_CI.Properties.RowNames = resExp.Properties.RowNames(idxexp);

end