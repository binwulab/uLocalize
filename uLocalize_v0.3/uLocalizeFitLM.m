function [final_pts, p, detection_res, fgImg]=uLocalizeFitLM(img, p, rough_pts, verbose)
%Objective: refine the rough points given in the pre-detection, one can
%also call this program to calculate intensity of given points 
% Currently supported method: Gaussian Mask or Gaussian Fit algorithm
%   In the future, can add other method p.method, for example
%   GaussianAnisotroicFit
%Input: 
%       img: a single 2D or 3D raw image
%       p: the parameter structure, which should include
%          {'width','thickness','sigma_xy','sigma_z','max_count',tol,'method'}
%          etc
%       rough_pts: the candidate points to fit 
%           [yc, xc, int] for 2D
%           [yc, xc, zc, int] for 3D
%       
%Output:
%   final_pts: 
%       for GaussianMask2D: [yc,xc,N0,Chi2]
%       for GaussianFit2D: [yc,xc,N0,sigma_xy]
%       for GaussianMask3D: [yc,xc,zc,N0,err]
%       for GaussianFit3D: [yc,xc,zc,N0,sigma_xy,sigma_z]
%   fgImg: the background image generated during the fit, could be useful for the TS quantification
%
%History   
%   B.W., Jan, 2021
%   B.W., Feb, 2021: change bgImg to fgImg

npts=size(rough_pts,1);
if npts == 0
    final_pts = [];
    return;
end
if nargin<4
    verbose=1;
end
%ROISizeOption='small';
%generate the foreground image
fgImg=genFGImageAroundPoints(img,rough_pts,p.cutwidth);    %Note: I removed all candidate spots from the image here, which will be used in the local background correction
if verbose
    disp('uLocalizeFit: beginning ...'); 
    tic;
end
if p.numdim==2
    sigma=p.sigma_xy;
else
    sigma=[p.sigma_xy, p.sigma_z];
end
if strcmpi(p.method, 'GaussianFit')
    final_pts=zeros(npts,p.numdim*2); %[y,x,z,I,sigma_xy,sigma_z]
else
    final_pts=zeros(npts,p.numdim+2); %[y,x,z,I,err0]
end
%% fit the points one by one
final_pts = uLocalizeFitSingleWrapper(img, fgImg, p, rough_pts);
valid_pts=find(~isnan(final_pts(:,1)));
final_pts=final_pts(valid_pts,:);

%% clean up array
if p.numdim == 2
    [ny,nx]=size(img);
    [~,final_pts,n_wrong, n_double]=clean_up_spots_array_2D_clean3(final_pts,ones(1,2),ny,nx,p.ROIsize);
else
    [ny,nx,nz]=size(img);
    [~,final_pts,n_wrong, n_double]=clean_up_spots_array_3D_clean3(final_pts,ones(1,3),ny,nx,nz,p.ROIsize);
end
if verbose 
    disp(['uLocalizeFitLM: eliminated ' num2str(n_wrong) ' wrong / ' num2str(n_double) ' double identifications. there remain ',num2str(size(final_pts,1)), 'spots on image.  Time spent: ' num2str(toc), ' seconds']);
end
detection_res=struct('rough_pts', rough_pts, 'final_pts', final_pts);
end