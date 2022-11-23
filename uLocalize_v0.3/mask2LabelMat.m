function [L, nL]=mask2LabelMat(mask, maskType)
%Obj: calculate the label matrix from a mask matrix
%Input
%   mask: a 2D mask matrix, each cell is a different value
%   maskType: 0: use the std of the whole image; 
%             1: each cell uses its own std; Default
%             2: all cells use the same std

sz=size(mask);
nY=sz(1);
nX=sz(2);
if maskType==0  %The whole image is used
    L=ones(nY, nX, 'uint16');
    nL=1;
elseif maskType==1  %each cell is calculated separately
    maxMask=max(mask(:));
    if maxMask > 1    %The mask should be Label matrix already
        L=uint16(mask);
        nL=maxMask;
    else
        [L,nL]=bwlabel(mask,4); %4-connected objects found in the mask
        L=uint16(L);
        if nL == 0  %mask is 0 everywhere
            disp('calcImgStd: empty mask');
            L=[];
            return
        end
    end
else    %All cells are used at the same time
    L=uint16(mask);
    L(L>0)=1;  %In this case, converting the label matrix to a binary mask
    nL=1;   
end
