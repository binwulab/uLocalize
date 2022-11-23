function cc=predetectCC(stack, thresh, maskType, mask,  areaRange, eccRange, verbose)
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
%   areaRange: [amin, amax], Area range
%   eccRange: [emin, emax], eccentricity range
%   ROISize: ROI size to determine whether it is the local maiximum
%Output
%   cc: the connected component structure, besides, NumObjects, ImageSize, Connectivity, PixelIdxList
%       nObjects: [nSingle, nMultiple, nZeros]
%       ObjCtr: the center of each object
%History   
%   B.W., Jan, 2021

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
%% set up the default value of input parameters
narg=nargin; 
if narg < 2 
    disp('predetectLM(stack, thresh)');
end
if narg<3;   maskType=1; end
if narg<4;   mask=ones(nY, nX); end

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

amin=areaRange(1);
amax=areaRange(2);
emin=eccRange(1);
emax=eccRange(2);
e=get_objects_planar_excentricity(cc);
% f0=a<amin & a==1;
% fb=a<amin & a>1;
% f1=a>=amin & a<=amax & e>=emin & e<=emax;
% fm=true(size(f1)) & (~f1) & (~fb) & (~f0);
%rearrange the objects to single, multiple and zeros
f0 = (a==1);
f1= a>=amin & a<=amax & e>=emin & e<=emax;
fm = a> amax & e>emax;
fb= true(size(f1)) & (~f1) & (~fm) & (~f0);
indfb=find(fb);
indf1=find(f1);
indf2=find(fm);
cc.PixelIdxList=cc.PixelIdxList([indf1,indf2, indfb]);
nObjects=[numel(indf1), numel(indf2), numel(indfb)];
cc.nObjects=nObjects;
cc.NumObjects=sum(cc.nObjects);
if verbose
    disp(['detected: ', num2str(cc.nObjects(1)), ' singles; ', num2str(cc.nObjects(2)), ...
        ' multiples; ', num2str(cc.nObjects(3)), ' bg spot; and removed ', num2str(numel(find(f0))), ...
        ' spots; maximum area: ', num2str(max(a))]);
    subplot(2,1,1);
    hist(a, 50); title('area')
    subplot(2,1,2);
    hist(e,50); title('eccentricity');
end

ObjCtr = regionprops(cc,'Centroid');
ObjCtr = reshape([ObjCtr.Centroid],nd, cc.NumObjects)';
ObjCtr(:,1:2) = [ObjCtr(:,2),ObjCtr(:,1)];
cc.ObjCtr=ObjCtr;
cc.threshInt=threshInt;
end