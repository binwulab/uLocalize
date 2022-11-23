function If = bpass_filter_2D_fast_memory_safe2(I,cutxy,dx,sigma_xy)
%performns a 2d (xy) bandpass filter on the 2d image I and stores it in Ismooth.

%input parameters: 
%dx is the voxel dimensions in nm (dx= dy assumed).
%sigma_xy is the PSF width in nm (assumed gaussian).
%cutxy is the lower cutoff frequency (fmin) of the bandpass filter in the xy dimensions, in units of the psf width.
%the higher cutoff frequency fmax is taken equal to 1.

tic;
[nx,ny] = size(I);
%removing the homogenous background 

if nx*ny > 696*520
    %I divide the large images into four quadrants to save memory
    npts = 11; % size of the kernel used in bpass_filter_2D_fast_on3Dimage; 
    %used to make enlarge each quadrant so that artifactual border regions
    %can then be excluded from the final image
    
    np2 = floor(npts/2);
    nx2 = round(nx/2); ny2 = round(ny/2);
    
    blacklevel = median_stack(I);
    
    I1 = I(1:nx2+np2,1:ny2+np2);
    I1 = I1 - blacklevel;
    I1 = bpass_filter_2D_fast_on2Dimage2(I1,cutxy,dx,sigma_xy,'same with zeros');
    If(1:nx2,1:ny2) = I1(1:nx2,1:ny2);
    clear('I1');
    
    I2 = I(nx2-np2+1:nx,1:ny2+np2);
    I2 = I2 - blacklevel;
    I2 = bpass_filter_2D_fast_on2Dimage2(I2,cutxy,dx,sigma_xy,'same with zeros');
    If(nx2+1:nx,1:ny2) = I2(np2+1:nx2+np2,1:ny2);
    clear('I2');
    
    I3 = I(1:nx2+np2,ny2-np2+1:ny);
    I3 = I3 - blacklevel;
    I3 = bpass_filter_2D_fast_on2Dimage2(I3,cutxy,dx,sigma_xy,'same with zeros');
    If(1:nx2,ny2+1:ny) = I3(1:nx2,np2+1:ny2+np2);
    clear('I3');
    
    I4 = I(nx2-np2+1:nx,ny2-np2+1:ny);
    I4 = I4 - blacklevel;
    I4 = bpass_filter_2D_fast_on2Dimage2(I4,cutxy,dx,sigma_xy,'same with zeros');
    If(nx2+1:nx,ny2+1:ny) = I4(np2+1:nx2+np2,np2+1:ny2+np2);
    clear('I4');  
    
else
    If = bpass_filter_2D_fast_on2Dimage2(I,cutxy,dx,sigma_xy,'same with zeros');
end
t = toc;
disp(['smoothing complete after ' num2str(t) ' s']);

clear('npts','nx','ny','nz','blacklevel','npts',...
    'np2','nx2','ny2','t');
end

function [If hxy] = bpass_filter_2D_fast_on2Dimage2(I,cutxy,dx,sigma_xy,shape)
%performns a 2d (xy) bandpass filter on the 3d image I and stores it in Ismooth.
%3d filter is too computationnally intensive and doesnt yield better image
%detection

%input parameters: 
%dx dy dz are the voxel dimensions in nm (dy not used, dx
%= dy assumed at this point).
%sigma_xy and sigma_z are the dimensions in nm of the PSF width.
%cutxy and cutz are the lower cutoff frequency (fmin) in the xy and z
%dimensions resp., in units of the psf width.
%the higher cutoff frequency fmax is taken equal to 1.

%outputs in addition of the bpassed image the filters in the frequency
%domain in the lateral plane 

npts = 11;  %size of the kernel; should be an odd integer; the largest the more computationnally intensive
[nx,ny,nz] = size(I);
%%%%%%% designing a 2D filter kernel for xy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmin = dx/(cutxy*sigma_xy);
fmax = 0.5;
[f1,f2] = freqspace(npts,'meshgrid');
%f1 and f2 are arrays of size npts, centered around zero, with values that
%are multiples of 2/npts.
Hd = ones(npts); 
r = sqrt(f1.^2 + f2.^2);
Hd((r<fmin)|(r>fmax)) = 0;

%plot the desired frequency response
%colormap(jet(64))
%figure;
%mesh(f1,f2,Hd);

%% Design the 2D filter h that approximate the desired freq response  
%using a 1D Hamming window.
hxy = fwind1(Hd,hamming(npts));
clear('Hd');
%plot the actual frequency response of the filter
%figure;
%freqz2(hxy);

%% abandoned
%get the desired frequency reponse in z:
%fmin = dz/(cutz*sigma_z);
%[f1,f2] = freqspace(npts,'meshgrid');
%f1 and f2 are arrays of size npts, centered around zero, with values that
%are multiples of 2/npts.
%Hz = ones(npts); 
%r = sqrt(f1.^2 + f2.^2);
%Hz((r<fmin)|(r>fmax)) = 0;

%% convolution in xy
if strcmp(shape,'same with zeros')
    If = zeros(nx,ny);
    If(ceil(npts/2):nx-floor(npts/2),ceil(npts/2):ny-floor(npts/2)) = convn(double(I),hxy,'valid');
elseif strcmp(shape,'valid')
    If(ceil(npts/2):nx-floor(npts/2),ceil(npts/2):ny-floor(npts/2)) = convn(double(I),hxy,'valid');
end
%[nx ny nz] = size(I);
%[newnx, newny, newnz] = size(Iconv);
%Ismooth = zeros(nx,ny,nz);
%diffx = (nx - newnx)/2; diffy = (ny - newny)/2; diffz = (nz - newnz)/2;

%Ismooth(diffx+1:diffx+newnx,diffy+1:diffy+newny,diffz+1:diffz+newnz) = Iconv(1:newnx,1:newny,1:newnz);

clear('nx','ny','nz','npts','fmin','fmax','f1','f2','Hd','r','Hz',...
    'diffx','diffy','diffz','hxyz');
end