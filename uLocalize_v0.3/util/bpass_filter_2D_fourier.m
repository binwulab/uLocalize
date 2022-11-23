function [If bp_filter] = bpass_filter_2D_fourier(I,cutxy,sigma_xy,l_hi,fw)


%performns a fourier space 2d (xy) bandpass filter on an image I and stores it into If.
%The filter shape in fourier space is stored in the image bp_filter

%cutoff frequencies:
%l_hi: lengthscale corresponding to the high frequency cutoff in pixel units.

%the lengthscale corresponding to the low frequency cutoff is defined as:
% cutxy*sigma_xy     (for single spot detection, sigma_xy is the PSF
% width in pixels), so cutxy is the cutoff length in PSF width
% units

%fw is the steepness of the filter in frequency space in pix^-1

%filt_order is the order of the filter (the higher the order, the steeper
%the filter
%the filter shape is a fermi dirac distribution



nx = size(I,1);
ny = size(I,2);

%fourier transform and centering of the spectrum
fI = fft2(I,nx,ny);
fI = fftshift(fI); % this puts the 0th order in the center of the image

%the array of the frequencies corresponding to the datapoints of fI
fx = 1/double(nx)*( (1:nx)  - 1 - floor(nx/2) );
fy = 1/double(ny)*( (1:ny)  - 1 - floor(ny/2) );

[xfx yfy] = meshgrid(fy, fx);
clear('fx','fy');

%design the filter:
% cutoff lengths in real space (in pix units)
l_lo = cutxy*sigma_xy;

%conversion in frequency space
fc_hi = 1/l_hi;
fc_lo = 1/l_lo;


r = sqrt(xfx.^2 + yfy.^2);
clear('xfx','yfy');
bp_filter = 1./(  (1+exp((r-fc_hi)./fw)) .* (1+exp((fc_lo - r)./fw))     ) ;
clear('r');

%apply the filter to the fourier transform of the image:
fIbp = fI .* bp_filter;
clear('fI');

%reconstruct the image in real space
fIbp = ifftshift(fIbp);
If = ifft2(fIbp);

clear('fIbp');


end