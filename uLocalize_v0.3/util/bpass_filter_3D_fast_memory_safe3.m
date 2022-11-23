function If = bpass_filter_3D_fast_memory_safe3(I,cutxy,sigma_xy,sigma_z,l_hi,filter_width,numdim)

%performns a fourier-based 3d (xy) bandpass filter on the 3d image I and stores it in Ismooth.

%input parameters: 
%dx dz are the voxel dimensions in nm (dy not used, dx= dy assumed).
%sigma_xy is the PSF width in nm.
%cutxy is the lower cutoff frequency (fmin) of the bandpass filter in the xy dimensions, in units of the psf width.
%the higher cutoff frequency fmax is taken equal to 1.
%l_hi,filter_width are filter parameters. Typically you can
%use l_hi = 3 pix, filter_width = 0.1 pix ^-1

tic;
[nx,ny,nz] = size(I);

if numdim ==2
    method = '2D';
elseif numdim == 3
    method = '3D';
end

if nx*ny > 512*512
    %I divide the large images into four quadrants to save memory
    nx2 = round(nx/2); ny2 = round(ny/2);
    if strcmp(method,'3D')
        %quadrant 1
        If(1:nx2,1:ny2,:)  = bpass_filter_3D_fourier_iso(I(1:nx2,1:ny2,:),cutxy,sigma_xy,sigma_z,l_hi,filter_width);

        %quadrant 2
        If(nx2+1:nx,1:ny2,:)  = bpass_filter_3D_fourier_iso(I(nx2+1:nx,1:ny2,:),cutxy,sigma_xy,sigma_z,l_hi,filter_width);

        %quadrant 3
        If(1:nx2,ny2+1:ny,:)  = bpass_filter_3D_fourier_iso(I(1:nx2,ny2+1:ny,:),cutxy,sigma_xy,sigma_z,l_hi,filter_width);

        %quadrant 4
        If(nx2+1:nx,ny2+1:ny,:)  = bpass_filter_3D_fourier_iso(I(nx2+1:nx,ny2+1:ny,:),cutxy,sigma_xy,sigma_z,l_hi,filter_width);
    elseif strcmp(method,'2D')
        If = zeros(nx,ny,nz);
        for i=1:nz
            If(:,:,i) = bpass_filter_2D_fourier(I(:,:,i),cutxy,sigma_xy,l_hi,filter_width);
        end
    end
else
    If  = bpass_filter_3D_fourier_iso(I,cutxy,sigma_xy,sigma_z,l_hi,filter_width);
end
t = toc;
disp(['smoothing complete after ' num2str(t) ' s']);

clear('npts','nx','ny','nz','blacklevel','npts',...
    'np2','nx2','ny2','t');

end






