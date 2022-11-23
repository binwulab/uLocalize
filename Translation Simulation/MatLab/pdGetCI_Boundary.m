function gfit = pdGetCI_Boundary(arginCI,N,IN,PD,PDnRMSE,arrayExp,arrayL,resExp,resGFit)

    gfit = resGFit;
    
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

    idx = find(nrmse(i1,:,:)-var<=arginCI*length(idxexp));
    [i2,i3] = ind2sub([N.koff,N.kon],idx);
    bd = boundary(i3,i2,0.1);
    z = zeros(length(bd),length(idxexp));

    for i = 1:length(idxexp)
        for j = 1:length(bd)
            z(j,i) = PD(i1,i2(bd(j)),i3(bd(j)),i,1);
        end
    end

    Boundary.ORF.IBID_MB = array2table([i3(bd),i2(bd),z],'VariableNames',{'X_ID','Y_IB','Z_MB_ST-AID','Z_MB_ORF1','Z_MB_ORF2','Z_MB_ORF3'});

    gfit.Boundary = Boundary;

end

