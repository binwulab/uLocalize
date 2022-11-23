function nrmse = pdGetRMSE(N,PD,arrayExp,arrayL)    

    nrmse = zeros(N.ki,N.koff,N.kon,size(arrayExp,1),4);
    
    for i = 1:size(arrayExp,1)       
        
        for j = 1:3
            wgt = arrayExp(i,2*(j-1)+2);
            nrmse(:,:,:,i,j) = (sqrt((PD(:,:,:,arrayL(i),2*(j-1)+1)-arrayExp(i,2*(j-1)+1)).^2))./wgt;
        end

        nrmse(:,:,:,i,4) = nrmse(:,:,:,i,1)+nrmse(:,:,:,i,2)+nrmse(:,:,:,i,3);

    end
    
end