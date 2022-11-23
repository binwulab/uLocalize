function cc=predetectGranule(stack, thresh, maskType, mask, amin, emax, verbose)
%Obj: given an image of 2D or 3D, and a threshold level in std, find object
%above a threshold, which is determined either for each cell or all cells,
%using connected commponent method
%Input
%   stack: raw image
%   threshLevelInSTD: the threshold in std
%   mask: a 2D label matrix, each cell is a different value
%   maskType: 0: use the std of the whole image; 
%             1: each cell uses its own std; Default
%             2: all cells use the same std
%   amin: minimum area
%   emax, eccentricity max
%Output
%   cc: the connected component structure, 
%       ObjCtr: the center of each object
%History   
%   B.W., Mar, 2021

sz=size(stack);
if length(sz)==2
    nd=2;
    nZ=1;
elseif length(sz) == 3
    nd=3;
    nZ=sz(3);
else
    cc=[];
    disp('predetectWithThreshInSTD: incorrect size of image');
    return
end
nY=sz(1);
nX=sz(2);

[L, nL]=mask2LabelMat(mask, maskType);
clear('mask');
threshInt=calcImgThresh(stack, thresh, L, nL);

%% thresholding
bw=false(size(stack));
for l=1:nL %
    maskTmp=false(nY,nX);
    maskTmp(L==l)=true;
    bwTmp=repmat(maskTmp,1,1, nZ);
    bw=bw | (bwTmp & stack>threshInt(l));
end

%% detect objects using the thresholded image
if nd==2
    cc=bwconncomp(bw, 4); %for 2d, only consider immediate neighbor
else
    cc=bwconncomp(bw,18); %for 3d, consider next nearest neighbor, sqrt(2) away
end
%% classifiy the objects based on its properties: size and eccentricity
a=regionprops(cc,'area');
a=[a.Area];
e=get_objects_planar_excentricity(cc);
f1= a>amin & e<emax; 
cc.PixelIdxList=cc.PixelIdxList(f1);
if verbose
    disp(['detected: ', num2str(numel(cc.PixelIdxList)), ' granules; and removed ', num2str(cc.NumObjects-numel(cc.PixelIdxList))]);
end
cc.NumObjects=numel(cc.PixelIdxList);

ObjCtr = regionprops(cc,'Centroid');
ObjCtr = reshape([ObjCtr.Centroid],nd, cc.NumObjects)';
ObjCtr(:,1:2) = [ObjCtr(:,2),ObjCtr(:,1)];  %change the order of x,y
cc.ObjCtr=ObjCtr;
cc.threshInt=threshInt;
end