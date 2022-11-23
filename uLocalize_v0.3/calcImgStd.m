function [sd]=calcImgStd(stack, L, nL)
%Obj: calculate the std of a given 2d or 3d stack, with option to calculate
%std for individual cells
%Input
%   stack: raw image
%   L: the label matrix
%   nL: the number of level in the label matrix. returned so no need to
%       calculate again
%Output
%   sd: for mask=0 or 2, return one sd value, for mask=1, return sd for
%       each cell

sz=size(stack);
if nargin<3 
    nL=max(L(:));
end
if length(sz)==2
    nd=2;
    nZ=1;
elseif length(sz) == 3
    nd=3;
    nZ=sz(3);
else
    disp('calcImgStd: can only handle 2d or 3d images');
    sd=[];
    return
end
nY=sz(1);
nX=sz(2);

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

end