function threshInt=calcImgThresh(stack, thresh, L, nL)
%Obj: calculate the std of a given 2d or 3d stack, with option to calculate
%std for individual cells
%Input
%   stack: raw image
%   thresh: the threshold structure, containing 'unit', 'level' field
%   L: the label matrix
%   nL: the number of level in the label matrix. returned so no need to
%       calculate again
%Output
%   sd: for mask=0 or 2, return one sd value, for mask=1, return sd for
%       each cell

sz=size(stack);
if nargin<2 || ~isstruct(thresh) || ~isfield(thresh, 'unit') || ~isfield(thresh,'level')
    disp('calcImgThresh: need threshold, dear');
    threshInt=[];
end
if nargin<3
    L=ones(sz(1:2), 'uint16');
end
if nargin<4 
    nL=max(L(:));
end
if length(sz)==2
    nZ=1;
elseif length(sz) == 3
    nZ=sz(3);
else
    disp('calcImgThresh: can only handle 2d or 3d images');
    threshInt=[];
    return
end

switch lower(thresh.unit)
    case 'sd'
        sd=zeros(1, nL);    %calculate std for each cell
        for l=1:nL %
            indx=find(L == l);  %all points for the current ROI
            %Find the std for the current ROI
            m1=zeros(1,nZ);
            m2=zeros(1,nZ);
            for i=1:nZ
                arr=stack(:,:,i);
                m1(i)=mean(arr(indx));
            end
            m1=mean(m1);
            for i=1:nZ
                arr=stack(:,:,i);
                m2(i)=mean((arr(indx)-m1).^2);
            end
            sd(l)=sqrt(mean(m2));
        end
        threshInt=sd*thresh.level;
    case 'absolute'
        threshInt=repmat(thresh.level,1,nL);    %
    case 'automatic'    %to be done
    otherwise
        disp('calcImgThresh: unknown threshold unit');
        threshInt=[];
        return
end
end