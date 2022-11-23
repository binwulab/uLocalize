function If = bpass_filter_3D_fast_memory_safe2(I,cutxy,dx,sigma_xy,cutz,dz,sigma_z,method)
%performns a 2d (xy) bandpass filter on the 3d image I and stores it in Ismooth.
%3d filter is too computationnally intensive and doesnt yield better image
%detection
%input parameters: 
%dx dz are the voxel dimensions in nm (dy not used, dx= dy assumed).
%sigma_xy is the PSFG width in nm.
%cutxy is the lower cutoff frequency (fmin) of the bandpass filter in the xy dimensions, in units of the psf width.
%the higher cutoff frequency fmax is taken equal to 1.

tic;
[nx,ny,nz] = size(I);
npts = 12; % size of the kernel used in bpass_filter_2D_fast_on3Dimage; 
    %used to make enlarge each quadrant so that artifactual border regions
    %can then be excluded from the final image
np2 = floor(npts/2);    %make that number odd

if nx*ny > 696*520
    %I divide the large images into four quadrants to save memory
    
    
    
    nx2 = round(nx/2); ny2 = round(ny/2);

    if strcmp(method,'2D filter')
        If(1:nx2,1:ny2,:) = bpass_filter_2D_fast_on3Dimage2(I(1:nx2+2*np2,1:ny2+2*np2,:),cutxy,dx,sigma_xy,'quadrant 1',2*np2+1);
    else
        If(1:nx2,1:ny2,:) = bpass_filter_3D_gaussian_on3Dimage2(I(1:nx2+np2,1:ny2+np2,:),cutxy,dx,sigma_xy,cutz,dz,sigma_z,'same with zeros',2*np2+1);
    end
    m1 = median_stack(If(1:nx2,1:ny2,:));
    
    if strcmp(method,'2D filter')
        If(nx2+1:nx,1:ny2,:) = bpass_filter_2D_fast_on3Dimage2(I(nx2-2*np2+1:nx,1:ny2+2*np2,:),cutxy,dx,sigma_xy,'quadrant 2',2*np2+1);
    else
        If(nx2+1:nx,1:ny2,:) = bpass_filter_3D_gaussian_on3Dimage2(I(nx2-np2+1:nx,1:ny2+np2,:),cutxy,dx,sigma_xy,cutz,dz,sigma_z,'same with zeros',2*np2+1);
    end
    m2 = median_stack(If(nx2+1:nx,1:ny2,:));
    
    if strcmp(method,'2D filter')
        If(1:nx2,ny2+1:ny,:) = bpass_filter_2D_fast_on3Dimage2(I(1:nx2+2*np2,ny2-2*np2+1:ny,:),cutxy,dx,sigma_xy,'quadrant 3',2*np2+1);
    else
        If(1:nx2,ny2+1:ny,:) = bpass_filter_3D_gaussian_on3Dimage2(I(1:nx2+np2,ny2-np2+1:ny,:),cutxy,dx,sigma_xy,cutz,dz,sigma_z,'same with zeros',2*np2+1);
    end
    m3 = median_stack(If(1:nx2,ny2+1:ny,:));
    
    if strcmp(method,'2D filter')
        If(nx2+1:nx,ny2+1:ny,:) = bpass_filter_2D_fast_on3Dimage2(I(nx2-2*np2+1:nx,ny2-2*np2+1:ny,:),cutxy,dx,sigma_xy,'quadrant 4',2*np2+1);
    else
        If(nx2+1:nx,ny2+1:ny,:) = bpass_filter_3D_gaussian_on3Dimage2(I(nx2-np2+1:nx,ny2-np2+1:ny,:),cutxy,dx,sigma_xy,cutz,dz,sigma_z,'same with zeros',2*np2+1);
    end
    m4 = median_stack(If(nx2+1:nx,ny2+1:ny,:));
    clear('I4');  
    
    %If = If - mean([m1;m2;m3;m4]);
    clear('m1','m2','m3','m4');
else
    if strcmp(method,'2D filter')
        If = bpass_filter_2D_fast_on3Dimage2(I,cutxy,dx,sigma_xy,'same with zeros',2*np2+1);
    else
        If = bpass_filter_3D_gaussian_on3Dimage2(I,cutxy,dx,sigma_xy,cutz,dz,sigma_z,'same with zeros',2*np2+1);
    end
end
t = toc;
disp(['smoothing complete after ' num2str(t) ' s']);

clear('npts','nx','ny','nz','blacklevel','npts',...
    'np2','nx2','ny2','t');

end

function [If hxy] = bpass_filter_2D_fast_on3Dimage2(I,cutxy,dx,sigma_xy,shape,npts)
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

  %size of the kernel; should be an odd integer; the largest the more computationnally intensive
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
clear('Hd','fmin','fmax','f1','f2','Hd','r');
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
    If = zeros(nx,ny,nz);
    for i=1:nz
        If(ceil(npts/2):nx-floor(npts/2),ceil(npts/2):ny-floor(npts/2),i) = convn(double(I(:,:,i)),hxy,'valid');
    end

elseif strcmp(shape,'valid')
    for i=1:nz
        If(ceil(npts/2):nx-floor(npts/2),ceil(npts/2):ny-floor(npts/2),i) = convn(double(I(:,:,i)),hxy,'valid');
    end

elseif strcmp(shape,'quadrant 1')    
    If = zeros(nx-npts+1,ny-npts+1,nz);
    for i=1:nz
        Ic = convn(double(I(:,:,i)),hxy,'valid');
        If(floor(npts/2)+1 : nx - 2*floor(npts/2),floor(npts/2)+1 : ny - 2*floor(npts/2),i) = Ic(1:nx - 3*floor(npts/2),1:ny - 3*floor(npts/2));
    end
    
    
elseif strcmp(shape,'quadrant 2')    
    If = zeros(nx-npts+1,ny-npts+1,nz);
    for i=1:nz
        Ic = convn(double(I(:,:,i)),hxy,'valid');
        If(1 : nx - 3*floor(npts/2),floor(npts/2)+1 : ny - 2*floor(npts/2),i) = Ic(floor(npts/2)+1:nx - 2*floor(npts/2),1:ny - 3*floor(npts/2));
    end
    
    
elseif strcmp(shape,'quadrant 3')    
    If = zeros(nx-npts+1,ny-npts+1,nz);
    for i=1:nz
        Ic = convn(double(I(:,:,i)),hxy,'valid');  
        If(floor(npts/2)+1 : nx - 2*floor(npts/2),1 : ny - 3*floor(npts/2),i) = Ic(1:nx - 3*floor(npts/2),floor(npts/2)+1:ny - 2*floor(npts/2));
    end
    
elseif strcmp(shape,'quadrant 4')    
    If = zeros(nx-npts+1,ny-npts+1,nz);
    for i=1:nz
        Ic = convn(double(I(:,:,i)),hxy,'valid');  
        If(1 : nx - 3*floor(npts/2),1 : ny - 3*floor(npts/2),i) = Ic(floor(npts/2)+1:nx - 2*floor(npts/2),floor(npts/2)+1:ny - 2*floor(npts/2));
    end 
    
end
%[nx ny nz] = size(I);
%[newnx, newny, newnz] = size(Iconv);
%Ismooth = zeros(nx,ny,nz);
%diffx = (nx - newnx)/2; diffy = (ny - newny)/2; diffz = (nz - newnz)/2;

%Ismooth(diffx+1:diffx+newnx,diffy+1:diffy+newny,diffz+1:diffz+newnz) = Iconv(1:newnx,1:newny,1:newnz);

clear('nx','ny','Ic','nz','npts','fmin','fmax','f1','f2','Hd','r','Hz',...
    'diffx','diffy','diffz','hxyz');
end

function [If kernel]= bpass_filter_3D_gaussian_on3Dimage2(I,cutxy,dx,sigma_xy,cutz,dz,sigma_z,shape)
%performns a 3d (xyz) gaussian convolution filter on the 3d image I and stores it in
%If.

%input parameters: 
%dx dy dz are the voxel dimensions in nm (dy not used, dx
%= dy assumed at this point).
%sigma_xy and sigma_z are the dimensions in nm of the PSF width.
%cutxy and cutz are the lower cutoff frequency (fmin) in the xy and z
%dimensions resp., in units of the psf width.
%the higher cutoff frequency fmax is taken equal to 1.

%outputs in addition of the bpassed image the convolution kernel
kernel = gaussian_kernel2(cutxy*sigma_xy/dx,cutz*sigma_z/dz,'peak');
If = imfilter(I,kernel,'replicate');


clear('nx','ny','nz','npts','fmin','fmax','f1','f2','Hd','r','Hz',...
    'diffx','diffy','diffz','hxyz');
end

function kernel = gaussian_kernel2(sigma_xieta,sigma_zeta,normalization_method)

nz = ceil(3*sigma_zeta); % at x = 3 sigma, exp(-x*x/(2*sigma*sigma)) = exp(-4.5) = 0.0111
nxy = ceil(3*sigma_xieta);

sigma2_xieta = sigma_xieta*sigma_xieta;
sigma2_zeta = sigma_zeta*sigma_zeta;
aux = zeros(2*nxy+1,2*nxy+1,2*nz+1);

for i = -nxy:nxy
    for j = -nxy:nxy
        for k = -nz:nz
            r2 = i.*i + j.*j;
            z2 = k.*k;
            aux(i+nxy+1,j+nxy+1,k+nz+1) = exp( -0.5 * (r2/sigma2_xieta + z2/sigma2_zeta) );
        end
    end
end

% normalize in the manner specified by the input argument 'normalization_method'
switch normalization_method
    
    case 'peak'
        aux0 = max(aux(:));
        
    case 'energy'
        aux0 = sum(aux(:));
        
    otherwise
        error(['normalization_method = ',normalization_method,' is not a valid option in gaussian_kernel !']);
        
end

kernel = aux/aux0;
clear('nz','nxy','sigma2_xieta','sigma2_zeta','aux','i','j','k','r2','z2','aux0');
end




