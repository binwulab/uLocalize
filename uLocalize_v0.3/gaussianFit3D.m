function [pt, status]=gaussianFit3D(Im,pt0, fitType, bw)
%The Gaussian Fit Wrapper program for cropped image
%   Imax*exp(-((x-xc)^2+(y-yc)^2)/(2*sigma^2))+bg
%Input
%       Im: the cropped image
%       pt0: the candidate point
%       sigma_xy: the psf sigma
%       fitType: 'standard gaussian'
%                'integrated gaussian'
%                'integrated gaussian no bg';
%       bw: uint8, the mask of other obj: >=1 if the pixel belongs to another object,it will
%           be excluded. The current obj, however, should be set 0. 

%Output
%       pt: the output [yc, xc, zc, Int, sxy, sz, bg]
%
%History
%   BW, Jan, 2021
%   BW, Feb, 2021: Add parameter bw to take into account overlapping objects
sz=size(Im);
if nargin<3 %Default algorithm
    fitType='integrated gaussian'; 
end 
if nargin<4 
    bw=zeros(sz,'uint8'); 
end
yc=ceil(pt0(1));
xc=ceil(pt0(2));
zc=ceil(pt0(3));
ny=sz(1);
nx=sz(2);
nz=sz(3);
data=double(Im);
[X,Y,Z]=meshgrid(1:nx, 1:ny, 1:nz);
indx=find(bw<1);
status=0;
if numel(indx)<8
    disp('gaussianFit3D: insufficient number of data points for fitting');
    status=-1;
    pt=-1+zeros(1,7);
    return;
end
pts=[Y(indx),X(indx),Z(indx)];
data=data(indx);
bg=min(data(:));
Imax=max(data(:))-bg;

%% setting up the parameters for the fitter
sxy=1;
sz=1;
Coeffs = [  Imax,   bg,         sxy,   sz,      yc,       xc,       zc];
lb =     [  0,      -Imax,      0.1,   0.1,     0,        0,        0];
ub =     [3000*Imax, Imax+bg,   20,    50,      ny,      nx,        nz];
options = optimset('TolX',0.0001,'MaxIter',1000,'Display','off');

%% Calling the fitter
if strcmp(fitType,'integrated gaussian')
    [Gaussout, resnorm]=lsqcurvefit(@gauss_integrated3D, Coeffs, pts,data,lb,ub,options);
elseif strcmp(fitType,'standard gaussian')
    [Gaussout, resnorm]=lsqcurvefit(@gauss3D, Coeffs, pts,data,lb,ub,options);
elseif strcmp(fitType,'integrated gaussian no bg')
    [Gaussout, resnorm]=lsqcurvefit(@gauss_integrated3D_no_bg, Coeffs, pts,data,lb,ub,options);
end
%% Reassign the parameter in the order: yc, xc, zc, intensity, sxy, sz, bg
pt=[Gaussout(5), Gaussout(6), Gaussout(7), Gaussout(1), Gaussout(3), Gaussout(4), Gaussout(2)];
end


%% the standard gaussian
function I = gauss3D(Coeffs, pts)
%Here the Gaussian is not normalized. 
% gridSize = size(pts,2);
% r = pts(:,1:gridSize/3,:);
% c = pts(:,gridSize/3+1:2*gridSize/3,:);
% h = pts(:,2*gridSize/3+1:gridSize,:);
r=pts(:,1);
c=pts(:,2);
h=pts(:,3);
I = Coeffs(2) + Coeffs(1)*exp( -( (r-Coeffs(5)).^2 + (c-Coeffs(6)).^2 )/(2*Coeffs(3)^2) + (h-Coeffs(7)).^2/(2*Coeffs(4)^2) );
end



%% the gaussian that takes into account the integration of the intensity
% over each pixel
function I = gauss_integrated3D(Coeffs, pts)
%everything is in pixel units here
%origin of positions at the corner on the image 
%(i.e. centers of pixels have half integer vox_size values)
%Note: the CDF of Gaussian(mu,sigma) is (1+erf((x-mu)/sqrt(2)*sigma))/2. So
%the function is normalized and Imax here should be the total photons. 
Imax = Coeffs(1); bg = Coeffs(2); sxy=Coeffs(3); sz=Coeffs(4); xc=Coeffs(5); yc=Coeffs(6); zc=Coeffs(7);

% gridSize = size(pts,2);
% r = pts(:,1:gridSize/3,:);
% c = pts(:,gridSize/3+1:2*gridSize/3,:);
% h = pts(:,2*gridSize/3+1:gridSize,:);
r=pts(:,1);
c=pts(:,2);
h=pts(:,3);

sqrt2_sxy=sqrt(2)*sxy;
diffx1 =  ((double(r-1)) - xc)/sqrt2_sxy;
diffx2 =  (double(r) - xc)/sqrt2_sxy;
intensity1 = abs( erf( diffx1) - erf(diffx2) );

diffy1 =  ((double(c-1)) - yc)/sqrt2_sxy;
diffy2 =  (double(c) - yc)/sqrt2_sxy;
intensity2 = abs( erf( diffy1) - erf(diffy2) );

sqrt2_sz=sqrt(2)*sz;
diffz1 =  ((double(h-1)) - zc)/sqrt2_sz;
diffz2 =  (double(h) - zc)/sqrt2_sz;
intensity3 = abs( erf( diffz1) - erf(diffz2) ); 

intensity = Imax.*intensity1.*intensity2.*intensity3;
I = intensity / 8.0 + bg;
%Note: 
end



%% same as above with background fixed to zero
function I = gauss_integrated3D_no_bg(Coeffs, pts)
%everything is in pixel units here
%origin of positions at the corner on the image 
%(i.e. centers of pixels have half integer vox_size values)

Imax = Coeffs(1); %bg = Coeffs(2);
sxy=Coeffs(3); sz=Coeffs(4); xc=Coeffs(5); yc=Coeffs(6); zc=Coeffs(7);

% gridSize = size(pts,2);
% r = pts(:,1:gridSize/3,:);
% c = pts(:,gridSize/3+1:2*gridSize/3,:);
% h = pts(:,2*gridSize/3+1:gridSize,:);
r=pts(:,1);
c=pts(:,2);
h=pts(:,3);

sqrt2_sxy=sqrt(2)*sxy;
diffx1 =  ((double(r-1)) - xc)/sqrt2_sxy;
diffx2 =  (double(r) - xc)/sqrt2_sxy;
intensity1 = abs( erf( diffx1) - erf(diffx2) );

diffy1 =  ((double(c-1)) - yc)/sqrt2_sxy;
diffy2 =  (double(c) - yc)/sqrt2_sxy;
intensity2 = abs( erf( diffy1) - erf(diffy2) );

sqrt2_sz=sqrt(2)*sz;
diffz1 =  ((double(h-1)) - zc)/sqrt2_sz;
diffz2 =  (double(h) - zc)/sqrt2_sz;
intensity3 = abs( erf( diffz1) - erf(diffz2) ); 

intensity = Imax.*intensity1.*intensity2.*intensity3;
I = intensity / 8.0;

end