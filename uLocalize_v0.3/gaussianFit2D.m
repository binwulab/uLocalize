function [pt, status]=gaussianFit2D(Im,pt0, fitType,bw)
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
%       pt: the output [yc, xc, Int, sxy, bg]
%History
%   BW, Jan, 2021
%   BW, Feb, 2021: Add parameter bw to take into account overlapping objects

sz=size(Im);
if nargin<3 %default algorithm
    fitType='integrated gaussian';
end
if nargin<4 
    bw=zeros(sz,'uint8'); 
end
yc=ceil(pt0(1));
xc=ceil(pt0(2));

ny=sz(1);
nx=sz(2);
data=double(Im);
[X,Y]=meshgrid(1:nx, 1:ny);
indx=find(bw<1);
status=0;
if numel(indx)<6
    disp('gaussianFit2D: insufficient number of data points for fitting');
    status=-1;
    pt=-1+zeros(1,5);
    return;
end

pts=[Y(indx),X(indx)];
data=data(indx);
bg=min(data(:));
Imax=max(data(:))-bg;
if Imax == 0
    u1 =     max([max(data(:)),1]);
else
    u1 = 3000*(Imax);
end

if Imax + bg == 0
    u2 =     max(double(max(max(data))),1);
else
    u2 = Imax + bg;
end

sxy=1;
Coeffs = [  Imax,   bg,     sxy,   yc,        xc];
lb =     [  0,      -u2,    0.1,   0,         0];
ub =     [  u1,     u2,     10,    ny,        nx];


options = optimset('TolX',0.0001,'MaxIter',1000,'Display','off');
if strcmp(fitType,'standard gaussian')
    [Gaussout, resnorm]=lsqcurvefit(@gauss2D, Coeffs, pts,data,lb,ub,options);
elseif strcmp(fitType,'integrated gaussian')
    [Gaussout, resnorm]=lsqcurvefit(@gauss_integrated2D, Coeffs, pts,data,lb,ub,options);
elseif strcmp(fitType,'integrated gaussian no bg')
    [Gaussout, resnorm]=lsqcurvefit(@gauss_integrated2D_no_bg, Coeffs, pts,data,lb,ub,options);
end
pt=[Gaussout(4), Gaussout(5), Gaussout(1)*4, Gaussout(3), Gaussout(2)];
end
%% the standard gaussian
function I = gauss2D(Coeffs, pts)

% gridSize = size(pts,2);
% r = pts(: , 1 : gridSize/2 );
% c = pts(: , gridSize/2 + 1 :gridSize );

r=pts(:,1);
c=pts(:,2);

Imax = Coeffs(1); bg = Coeffs(2); sxy=Coeffs(3); yc=Coeffs(4); xc=Coeffs(5);

I = bg + Imax*exp( -( (r-yc).^2 + (c-xc).^2 )/(2*sxy^2)  );
end



%% the gaussian that takes into account the integration of the intensity
% over each pixel
function I = gauss_integrated2D(Coeffs, pts)
%everything is in pixel units here
%origin of positions at the corner on the image 
%(i.e. centers of pixels have half integer vox_size values)

Imax = Coeffs(1); bg = Coeffs(2); sxy=Coeffs(3); xc=Coeffs(4); yc=Coeffs(5);

% gridSize = size(pts,2);
% r = pts( : , 1 : gridSize/2 );
% c = pts( : , gridSize/2 + 1 : gridSize );
r=pts(:,1);
c=pts(:,2);

diffx1 =  (double(r-1)) - xc;
diffx1 = diffx1 ./ ( sqrt(2) * sxy );
diffx2 =  double(r) - xc;
diffx2 = diffx2 ./ ( sqrt(2) * sxy );

intensity1 = abs( erf( diffx1) - erf(diffx2) );

diffy1 =  (double(c-1)) - yc;
diffy1 = diffy1 ./ ( sqrt(2) * sxy );
diffy2 =  double(c) - yc;
diffy2 = diffy2 ./ ( sqrt(2) * sxy );

intensity2 = abs( erf( diffy1) - erf(diffy2) );

intensity = Imax.*intensity1.*intensity2;
I = intensity  + bg;

end




%% same as previous, but assumes bg is zero
function I = gauss_integrated2D_no_bg(Coeffs, pts) 
%everything is in pixel units here
%origin of positions at the corner on the image 
%(i.e. centers of pixels have half integer vox_size values)

Imax = Coeffs(1); %bg = Coeffs(2);
sxy=Coeffs(3); xc=Coeffs(4); yc=Coeffs(5); 

% gridSize = size(pts,2);
% r = pts( : , 1 : gridSize/2 );
% c = pts( : , gridSize/2 + 1 : gridSize );
r=pts(:,1);
c=pts(:,2);

diffx1 =  (double(r-1)) - xc;
diffx1 = diffx1 ./ ( sqrt(2) * sxy );
diffx2 =  double(r) - xc;
diffx2 = diffx2 ./ ( sqrt(2) * sxy );

intensity1 = abs( erf( diffx1) - erf(diffx2) );

diffy1 =  (double(c-1)) - yc;
diffy1 = diffy1 ./ ( sqrt(2) * sxy );
diffy2 =  double(c) - yc;
diffy2 = diffy2 ./ ( sqrt(2) * sxy );

intensity2 = abs( erf( diffy1) - erf(diffy2) );

intensity = Imax.*intensity1.*intensity2;
I = intensity;

end
