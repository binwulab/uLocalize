function [validPts, threshInt] = predetectLM(stack, thresh, maskType, mask,  ROISize, maxSpots, verbose)
%Obj: given an image of 2D or 3D, and a threshold level in std, find points
%above a threshold, which is determined either for each cell or all cells
%Input
%   stack: raw image
%   threshLevelInSTD: the threshold in std
%   mask: a 2D label matrix, each cell is a different value
%   maskType: 0: use the std of the whole image; 
%             1: each cell uses its own std; Default
%             2: all cells use the same std
%   ROISize: ROI size to determine whether it is the local maiximum
%Output
%   validPts: Valid local maximum points (Y,X,Z,Int), for 2D images, Z=1
%History   
%   B.W., Jan, 2021
%       Sep, 2021, Extract out the std calculation. allows absolute
%           threshLevel
sz=size(stack);
if length(sz)==2
    nd=2;
    nZ=1;
elseif length(sz) == 3
    nd=3;
    nZ=sz(3);
else
    validPts=zeros(0,4);
    disp('predetectWithThreshInSTD: incorrect size of image');
    return
end
nY=sz(1);
nX=sz(2);


%% set up the default value of input parameters
narg=nargin; 
if narg < 2 
    disp('predetectLM(stack, thresh)');
end
if narg<3;   maskType=1; end
if narg<4;   mask=ones(nY, nX); end
if narg<5;   ROISize=2; end
if narg<6;   maxSpots=50000; end

[L, nL]=mask2LabelMat(mask, maskType);
clear('mask');
threshInt=calcImgThresh(stack, thresh, L, nL);

pts=zeros(0,4);
for l=1:nL %
    indx=find(L == l);  %all points for the current ROI
    for i=1:nZ
        arr=stack(:,:,i);
        indxAboveThresh=indx(arr(indx)>=threshInt(l));
        [r,c]=ind2sub([nY,nX],indxAboveThresh);  %Note: points above thresh in the curr ROI
        np=length(r);
        if np>=1
            tmp=zeros(np,4);
            tmp(:,1)=r;
            tmp(:,2)=c;
            tmp(:,3)=i;
            tmp(:,4)=arr(indxAboveThresh);
            pts=[pts;tmp];  %Concatenate all detected points
        end
    end
end


%% remove redundant points
npts=size(pts,1);
validPts=zeros(npts,4);
nValid=1;
for i=1:npts    %go over each points to see whether it is the maximum around an ROI
    ymin=max(1,min(pts(i,1)-ROISize,nY));
    ymax=max(1, min(pts(i,1)+ROISize,nY));
    xmin=max(1,min(pts(i,2)-ROISize,nX));
    xmax=max(1, min(pts(i,2)+ROISize,nX));
    zmin=max(1,min(pts(i,3)-ROISize,nZ));
    zmax=max(1, min(pts(i,3)+ROISize,nZ));
    ROI=stack(ymin:ymax, xmin:xmax, zmin:zmax); %Note: we are using the original stack, not the masked one
    if max(ROI(:)) == pts(i,4)  %If the point is the local maximum, keep it
        validPts(nValid,:)=pts(i,:);
        nValid=nValid+1;
    end        
end
if nValid>1
    validPts=validPts(1:nValid-1,:);
else
    disp('Warning: predetectLM: no points detected');
    return;
end
%% Sort the spots according to the intensity
[~, indices]=sort(validPts(:,4), 'descend');    %sort pts in descending intensities
validPts=validPts(indices,:);
if size(validPts,1)>maxSpots    %If too many spots detected, only keeep the high intensity ones
    validPts=validPts(1:maxSpots,:);
    disp(['predetecLM: detected ', num2str(maxSpots), ' spots, maximum number achieved']);
elseif verbose
    disp(['predetecLM: detected ', num2str(size(validPts,1)), ' spots']);
end
% if nZ==1    %For 2D image, only return 3 columns [Y,X,Int]. 
%     validPts=validPts(:,[1,2,4]);
% end

%% building the arrays of coordinates and intensity of spots for subsequent fitting

end