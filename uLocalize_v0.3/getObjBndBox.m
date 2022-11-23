function [bndBox, pixelIdxList2, bw2]=getObjBndBox(pixelIdxList, imSize, ext, bw)
%objective: given an object in an image and the extension around the obj,
%calculate a bounding box, and the corresponding pixelIdxList in the
%cropped image, or the cropped mask matrix
%Input:
%   pixelIdxList: the linear index of the object
%   imSize: the size of the image [nY,nX,nZ]
%   ext: the number of pixel extended from the object, for 2D image, [ext_xy], for 3D, [ext_xy, ext_z] 
%   bw: the total mask matrix for original matrix, optional
%Output
%   bndBox: the bounding box around the object, extended from boundary of the object by [ext_xy, ext_z]
%   pixelIdxList2: the pixel list of the object in the new cropped image in the bounding box
%   bw2: the cropped mask image, or the generated cropped mask image
%History
%   BW, Jan 2021

numdim=numel(imSize);
if numdim == 3
    [yobj,xobj,zobj] = ind2sub(imSize,pixelIdxList);
    if numel(ext) == 1  %If only one extension para is given
        ext = [ext, ext];
    end
    %generating the bndBox for the object, leaving space 'ext' pix around it
    y1 = max([min(yobj)-ext(1),1]);
    y4 = min([max(yobj)+ext(1),imSize(1)]);
    x1 = max([min(xobj)-ext(1),1]);
    x4 = min([max(xobj)+ext(1),imSize(2)]);
    z1 = max([min(zobj)-ext(2),1]);
    z4 = min([max(zobj)+ext(2),imSize(3)]);
    bndBox=[y1,x1,z1; y4,x4,z4];    %The bounding box
    ny=y4-y1+1;
    nx=x4-x1+1;
    nz=z4-z1+1;
    if nargout > 1  %calculate only when necessary
        yobj=yobj-y1+1;     %compute the cropped object
        xobj=xobj-x1+1;
        zobj=zobj-z1+1;
        pixelIdxList2=(zobj-1)*ny*nx+(xobj-1)*ny+yobj;  %convert the pixel list to the coord in the cropped image
    end
    if nargout > 2          %generate the cropped mask file
        if nargin > 3       %if the input mask is given, just use it
            bw2=bw(y1:y4,x1:x4,z1:z4);
        else                %generate new mask file
            bw2=zeros(ny,nx,nz,'uint8');
            bw2(pixelIdxList2)=1;
        end
    end  
elseif numdim == 2
    [yobj,xobj] = ind2sub(imSize,pixelIdxList);
    %generating the bndBox for the object, leaving space ext around it
    y1 = max([min(yobj)-ext(1),1]);
    y4 = min([max(yobj)+ext(1),imSize(1)]);
    x1 = max([min(xobj)-ext(1),1]);
    x4 = min([max(xobj)+ext(1),imSize(2)]);
    bndBox=[y1,x1; y4,x4];    %The bounding box
    ny=y4-y1+1;
    nx=x4-x1+1;
    if nargout > 1  %calculate only when necessary
        yobj=yobj-y1+1; %compute the cropped object
        xobj=xobj-x1+1;
        pixelIdxList2=(xobj-1)*ny+yobj;  %convert the pixel list to the coord in the cropped image
    end
    %Note, even the 2nd parameter is not assigned, the nargout will be 3 if
    %the 3rd output arg is given. So pixelIdxList2 will be calculated here
    if nargout > 2      %generate the cropped mask file, 
        if nargin > 3   %if the input mask is given, just use it
            bw2=bw(y1:y4,x1:x4);
        else            %generate new mask file
            bw2=zeros(ny,nx,'uint8');
            bw2(pixelIdxList2)=1;
        end
    end  
end
end