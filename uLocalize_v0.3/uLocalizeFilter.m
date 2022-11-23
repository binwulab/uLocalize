function smooth=uLocalizeFilter(img, p)
%Objective: filter the image for pre-detection, right now, I only implement the
%bandpass filtering. Later will try LoG
if p.numdim == 2
    switch upper(p.filter.filterMethod)
        case 'BANDPASS'
            smooth=bpass_filter_2D_fourier(img,p.filter.nlo,p.sigma_xy,p.filter.nhi,p.filter.width);
        case 'LOGRAJ'
            smooth=filterLoGRaj(img, p.filter.filterSigma, p.filter.filterWindowSize);
        case 'LOGFFT'
            smooth=filterLoGFFT(img,p.filter.filterSigma);
            smooth=smooth-mean(smooth(:));
        otherwise
            disp('unknown filter in uLocalizeFilter');
            smooth=[];
            return
    end
else
    switch upper(p.filter.filterMethod)
        case 'BANDPASS'
            smooth=double(bpass_filter_3D_fast_memory_safe3(img,p.filter.nlo,p.sigma_xy,p.sigma_z,p.filter.nhi,p.filter.width,p.filter.numdim)); 
            smooth=smooth-mean(smooth(:));
        case 'LOGRAJ'
            smooth=filterLoGRaj(img, p.filter.filterSigma, p.filter.filterWindowSize);
        case 'LOGFFT'
            smooth=filterLoGFFT(img,p.filter.filterSigma);
            smooth=smooth-mean(smooth(:));
        otherwise
            disp('unknown filter in uLocalizeFilter');
            smooth=[];
            return
    end
end
end