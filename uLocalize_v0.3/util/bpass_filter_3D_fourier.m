function [If bp_filter] = bpass_filter_3D_fourier(I,cut,dx,dz,sigma_xy,sigma_z,l_hi,fw)


%performns a fourier space 2d (xy) bandpass filter on an image I and stores it into If.
%The filter shape in fourier space is stored in the image bp_filter

%cutoff frequencies:
%l_hi: lengthscale corresponding to the high frequency cutoff in pixel units.

%the lengthscale corresponding to the low frequency cutoff is defined as:
% cutxy*sigma_xy/dx     (for single spot detection, sigma_xy is the PSF
% width and dx the pixel size), so cutxy is the cutoff length in PSF width
% units

%fw is the steepness of the filter in frequency space in pix^-1


%the filter shape is the product of two fermi dirac distributions (one for
%the low pass, one for the high pass), to the power filter_order


[nx,ny,nz] = size(I);

%the array of the frequencies corresponding to the datapoints of the FT of
%the image
fx = 1/double(nx)*( (1:nx)  - 1 - floor(nx/2) );
fy = 1/double(ny)*( (1:ny)  - 1 - floor(ny/2) );

[xfx yfy] = meshgrid(fy, fx);
clear('fx','fy');

%design the filter:
% cutoff lengths in real space (in pix units)
l_lo_xy = cut*sigma_xy / dx;
l_lo_z = cut*sigma_z / dz;

%conversion in frequency space
fc_hi = 1/l_hi;
fc_lo_xy = 1/l_lo_xy;
fc_lo_z = 1/l_lo_z;

r = sqrt(xfx.^2 + yfy.^2 );
clear('xfx','yfy');
r = repmat(r,[1,1,nz]);

bp_filter = 1./(  (1+exp((r-fc_hi)./fw)) .* (1+exp((fc_lo_xy - r)./fw))   ) ;

fx = 1/double(nx)*( (1:nx)  - 1 - floor(nx/2) );
fy = 1/double(ny)*( (1:ny)  - 1 - floor(ny/2) );
fz = 1/double(nz)*( (1:nz)  - 1 - floor(nz/2) );

[~, ~, zfz] = meshgrid(fy, fx, fz);
r = abs(zfz);
clear('zfz');

bp_filter = bp_filter ./(  (1+exp((r-fc_hi)./fw)) .* (1+exp((fc_lo_z - r)./fw))   ) ;
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