function [pt,status]=gaussianMask2D(Im, pt0, sigma_xy, maxIter, tol, bw)
%function [pt,dist,it]=gaussianMask2D(Im, pt0, sigma_xy, maxIter, tol, bw)
%The Gaussain Mask Algorithm: given a cropped image containing a point
%source, use Gaussian Mask algorithm to find the center and intensity
%Input
%   Im: the cropped input image
%   pt0: the starting position [yc,xc,int,err0], the int is not used
%   cutwidth: the ROI size in pixel
%   maxIter: maximum iteration
%   tol: tolerance of the convergence
%   bw: uint8, the mask of other obj: >=1 if it belong to another object, will be excluded
%Output
%   pt: [yc, xc, int, err0], the only parameter that is used later
%   dist: the distance between successive iteration
%   it: the actual number of iteration
%History:
%   Apr 2021, change the output from [pt,dist,it] to [pt, status]

sz=size(Im);
if nargin<6
    bw=zeros(sz,'uint8');
end
ny=sz(1);
nx=sz(2);
Im=double(Im);
pt0=double(pt0);
y0=zeros(1,maxIter);
x0=zeros(1,maxIter);
N0=zeros(1,maxIter);

it=2;
y0(it)=pt0(1);
x0(it)=pt0(2);
y=pt0(1);
x=pt0(2);
dist=zeros(1,maxIter);

yp_min=1;
xp_min=1;
yp_max=ny;
xp_max=nx;
[xp,yp]=meshgrid(1:nx,1:ny);
% only consider the points that does not include other objects
indx=find(bw==0);
status=0;
if numel(indx)<3
    disp('gaussianMask2D: insufficient number of data points for fitting');
    status=-1;
    pt=-1+zeros(1,4);
    return;
end

xp=xp(indx);    %This also convert into linear indexing
yp=yp(indx);
Im=Im(indx);

tmp=1+tol;
sqrt2_sxy=sqrt(2)*sigma_xy; %pre-calculate some constant, so no need to re-calculate again and again
while (it<=maxIter && tmp>tol) %The actual Gaussian Mask Iteration 
    diffx1 = (xp-1-x)/sqrt2_sxy;
    diffx2 = (xp-x)/sqrt2_sxy;
    
    diffy1 = (yp-1-y)/sqrt2_sxy;
    diffy2 = (yp -y)/sqrt2_sxy;
    
    intensity = abs( erf( diffx1) - erf(diffx2) ) .* abs( erf( diffy1) - erf(diffy2) ); %Note: this is the integrated gaussian (which is error function) for each pixel, not directly the gaussian function
    
    intsum = intensity .* Im;
    sumsum = intensity.*intensity;
    sumx = (xp-0.5) .* intensity .* Im;
    sumy = (yp-0.5) .* intensity .* Im;

    intsum = abs( sum(intsum) );
    sumsum = abs( sum(sumsum) );
    sumx = sum(sumx);
    sumy = sum(sumy);

    if intsum <= 0 || sumsum == 0
        x0(it) = -1;
        y0(it) = -1;
        N0(it) = -1; 
    else
        x0(it) = double(sumx) / double(intsum);
        y0(it) = double(sumy) / double(intsum);
        N0(it) = 4*double(intsum) / double(sumsum); % this converts the previous value of N, Used for the integrated gaussian algorithm
        %(prefactor of a 2D erf function) into
        %the actual number of counts integrated over space
        if   (  ceil(x0(it)) > xp_max || ceil(x0(it)) < xp_min || ...
                ceil(y0(it)) > yp_max || ceil(y0(it)) < yp_min )
            x0(it) = -1;
            y0(it) = -1;
            N0(it) = -1;
        end
    end
   
    dist(it) = sqrt( (x-x0(it))^2 + (y-y0(it))^2);
    x = x0(it); y = y0(it);
    tmp = dist(it);
    if x0(1,it) == -1
        tmp = tol-1;
    end
    it = it+1;
end
if it >1
    if ceil(x0(it-1)) > xp_max || ceil(x0(it-1)) < xp_min ||  ceil(y0(it-1)) > yp_max || ceil(y0(it-1)) < yp_min 
        x0(1,it) = -1;
        y0(1,it) = -1;
        N0(1,it) = -1;
        dist(1,it) = -1;
        it = it+1;    
    end
    x0 = x0(1:it-1); y0 = y0(1:it-1);
    N0 = N0(1:it-1); dist = dist(1:it-1);
end

err0 = (N0(it-1)*intensity/4 - Im ).^2;
err0 = sum(err0(:))/numel(xp);% / ( (xp_max - xp_min+1) * (yp_max - yp_min+1) );
pt=[y0(end), x0(end), N0(end), err0];
end
