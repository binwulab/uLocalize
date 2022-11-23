function bgImg=genBGImgWrapper(par,Y,X,Z)
%Objective: Given the fitted parameter from background correction, generate
%the background image 
numdim=ndims(Y);
if numdim == 3 
    if nargin<4
        disp('GenBGCorrImgWrapper: 3D image, but only provided 2 coordinates, dear!');
        return;
    end
    switch numel(par)
        case 4  %local plane
            bgImg=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3)*(Z-0.5)+par(4);
        case 5  % local quad
            bgImg=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3)*(Z-0.5)+par(4)+par(5)*(Z-0.5).^2;
        case 1  % local median
            bgImg = zeros(size(Y))+par;
        otherwise %invalid
            disp('genBGImgWrapper 3D: : the number of parameters seems not right, dear!');
    end
else
    switch numel(par)
        case 3  %local plane
            bgImg=par(1)*(Y-0.5)+par(2)*(X-0.5)+par(3);    %bg correction in the imaging area through interpolation
        case 1  % local median
            bgImg = zeros(size(Y))+par;
        otherwise %invalid
            disp('genBGImgWrapper 2D: : the number of parameters seems not right, dear!');
    end
end
end