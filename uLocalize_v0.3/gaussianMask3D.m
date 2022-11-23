function [pt,status]=gaussianMask3D(Im,pt0, sigma, maxIter,tol,bw)
%function [pt,dist,it]=gaussianMask3D(Im,pt0, sigma, maxIter,tol,bw)
%The Gaussain Mask Algorithm: given a cropped image containing a point
%source, use Gaussian Mask algorithm to find the center and intensity
%Input
%   Im: the cropped input image
%   pt0: the starting position [yc,xc,zc,int,err0]. The int is not used
%   sigma: [sigma_xy, sigma_z]
%   maxIter: maximum iteration
%   tol: tolerance for the convergence
%   bw: uint8, the mask of other obj: >=1 if it belong to another object,will be excluded
%Output
%   pt: [yc, xc, zc, int, err], the only parameter that is used later
%   dist: the distance between successive iteration
%   it: the actual number of iteration
%History:
%   Apr 2021, change the output from [pt,dist,it] to [pt, status]

sz=size(Im);
if nargin<6
    bw=zeros(sz, 'unit8');
end
ny=sz(1);
nx=sz(2);
nz=sz(3);
pt0=double(pt0);
Im=double(Im);
y0=zeros(1,maxIter);
x0=zeros(1,maxIter);
z0=zeros(1,maxIter);
N0=zeros(1,maxIter);

it=1;
y0(it)=pt0(1);
x0(it)=pt0(2);
z0(it)=pt0(3);
N0(it)=0;

y=pt0(1);
x=pt0(2);
z=pt0(3);
dist=zeros(1,maxIter);

yp_min=1;
xp_min=1;
zp_min=1;
yp_max=ny;
xp_max=nx;
zp_max=nz;

[xp,yp,zp]=meshgrid(1:nx,1:ny,1:nz);
% only consider the points that does not include other objects
indx=find(bw==0);
status=0;
if numel(indx)<6
    disp('gaussianMask3D: insufficient number of data points for fitting');
    status=-1;
    pt=-1+zeros(1,5);
    return;
end

xp=xp(indx);    %This also convert into linear indexing
yp=yp(indx);
zp=zp(indx);
Im=Im(indx);

tmp=1+tol;
sqrt2_sxy=sqrt(2)*sigma(1); %pre-calculate some constant, so no need to re-calculate again and again
sqrt2_sz=sqrt(2)*sigma(2);
it=2;
while (it<=maxIter && tmp>tol)  %The actual Gaussian Mask Iteration 
    diffx1 = (xp-1-x)/sqrt2_sxy;
    diffx2 = (xp-x)/sqrt2_sxy;
    
    diffy1 = (yp-1-y)/sqrt2_sxy;
    diffy2 = (yp -y)/sqrt2_sxy;

    diffz1 = (zp-1-z)/sqrt2_sz;
    diffz2 = (zp -z)/sqrt2_sz;

    intensity = abs( erf( diffx1) - erf(diffx2) ).* abs( erf( diffy1) - erf(diffy2) ).* abs( erf( diffz1) - erf(diffz2) );
    
    intsum = intensity .* Im;
    sumsum = intensity.*intensity;
    sumx = (xp-0.5) .* intensity .* Im;
    sumy = (yp-0.5) .* intensity .* Im;
    sumz = (zp-0.5) .* intensity .* Im;

    intsum = abs(sum(intsum(:)));
    sumsum = abs(sum(sumsum(:)));
    sumx = sum(sumx(:));
    sumy = sum(sumy(:));
    sumz = sum(sumz(:));

    if intsum <= 0 || sumsum == 0
        x0(it) = -1;
        y0(it) = -1;
        z0(it) = -1;
        N0(it) = -1; 
    else
        x0(it) = double(sumx) / double(intsum);
        y0(it) = double(sumy) / double(intsum);
        z0(it) = double(sumz) / double(intsum);
        N0(it) = 8*double(intsum) / double(sumsum); % this converts the previous value of N, Used for the integrated gaussian algorithm 
        %(prefactor of a 2D erf function) into
        %the actual number of counts integrated over space
        if   (  ceil(x0(it)) > xp_max || ceil(x0(it)) < xp_min || ...
                ceil(y0(it)) > yp_max || ceil(y0(it)) < yp_min || ...
                ceil(z0(it)) > zp_max || ceil(z0(it)) < zp_min )
            x0(it) = -1;
            y0(it) = -1;
            z0(it) = -1;
            N0(it) = -1;
        end
    end
   
    dist(it) = sqrt( (x-x0(it))^2 + (y-y0(it))^2+(z-z0(it))^2);
    x = x0(it); y = y0(it); z=z0(it);
    tmp = dist(it);
    if x0(1,it) == -1
        tmp = tol-1;
    end
    it = it+1;
end

% if the spot converged outside of its ROI
if it >1
    if     (ceil(x0(it-1)) > xp_max || ceil(x0(it-1)) < xp_min ||  ...
            ceil(y0(it-1)) > yp_max || ceil(y0(it-1)) < yp_min || ...
            ceil(z0(it-1)) > zp_max || ceil(z0(it-1)) < zp_min )
        x0(it) = -1;
        y0(it) = -1;
        z0(it) = -1;
        N0(1,it) = -1;
        dist(1,it) = -1;
        it = it+1;    
    end
    x0 = x0(1:it-1); y0 = y0(1:it-1); z0=z0(1:it-1);
    N0 = N0(1:it-1); dist = dist(1:it-1);
end
err0=(N0(it-1)*double(intensity)/8.0 - Im).^2;
err0=sum(err0(:))/numel(xp);
pt=[y0(end), x0(end), z0(end), N0(end), err0];
end
