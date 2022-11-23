function [If bp_filter] = bpass_filter_2D_fourier_on_3D_stack(I,cut,sigma_xy,l_hi,fw)


%performns a fourier space 3d (xyz) bandpass filter on an image I and stores it into If.
%The filter shape in fourier space is stored in the image bp_filter

%cutoff frequencies:
%l_hi: lengthscale corresponding to the high frequency cutoff in pixel units.

%the lengthscale corresponding to the low frequency cutoff is defined as:
% cutxy*sigma_xy/dx     (for single spot detection, sigma_xy is the PSF
% width and dx the pixel size), so cutxy is the cutoff length in PSF width
% units

%fw is the steepness of the filter in frequency space in pix^-1

%filt_order is the order of the filter (the higher the order, the steeper
%the filter

%the filter has the shape of the shell of an almond.
%it is the product of two fermi dirac distributions (one for
%the low pass, one for the high pass)


%parameters that give good results for diffraction limited spots:
%cut = 3;

%sigma_xy = 2; (this should be the FWHM of the PSF in xy plane in pix units)
%sigma_z = 2; this should be the FWHM of the PSF in the axial direction in pix units)
%l_hi = 3:
%f_w = 0.1;

%size of the image
[nx,ny,nz] = size(I);

%the array of the frequencies corresponding to the datapoints of the FT of
%the image
fx = 1/double(nx)*( (1:nx)  - 1 - floor(nx/2) );
fy = 1/double(ny)*( (1:ny)  - 1 - floor(ny/2) );
fz = 1/double(nz)*( (1:nz)  - 1 - floor(nz/2) );

%collecting them in arrays that have the same size as the original stack
[r yfy zfz] = meshgrid(fy,fx,fz);
clear('zfz');
% cutoff lengths in real space (in pix units)
l_lo_xy = cut*sigma_xy;

%conversion in frequency space
fc_hi = 1/l_hi;
fc_lo_xy = 1/l_lo_xy;
clear('l_lo_xy','l_lo_z');

%building the low pass component of the filter
r = sqrt(r.^2 + yfy.^2); %radius in frequency space
clear('yfy');
bp_filter = 1./(1+exp((r - fc_hi)./fw)) ;
clear('rh');

%multiplying by the high pass component of the filter
[r yfy zfz] = meshgrid(fy,fx,fz);
clear('zfz');
r = sqrt( r.^2 + yfy.^2); %normalized 'isotropic' radius in frequency space
clear('yfy');
bp_filter = bp_filter ./((1+exp((fc_lo_xy-r)./fw))) ;
clear('r');

%fourier transform and centering of the spectrum
fI = fftn(I);
fI = fftshift(fI); % this puts the 0th order in the center of the image

%apply the filter to the fourier transform of the image:
fIbp = fI .* bp_filter;
clear('fI');

%reconstruct the image in real space
fIbp = ifftshift(fIbp);
If = ifftn(fIbp);
clear('fIbp');


end