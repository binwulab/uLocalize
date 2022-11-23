function fgImg=genFGImageAroundPoints(img, pts, cutwidth)
%Objective: given an image and a set of candidate points, generate a
%foreground image by adding 1 to the point and its surroundings in with the
%size of ROISize (by setting it to be nan). This fgImg can be used later for local background
%subtraction and fitting. This should work better if some points are close to each other. 
%Input: 
%   stack: 2D or 3D images
%   pts: the rough spots position determined from pre-detection [Y,X,Z,Int]
%   cutwidth: [cutwidth_xy, cutwidth_z], for 2D, only cutwidth_xy is provided
%Output
%   fgImg: the foreground image, each pixel represents the number of object it belongs to

sz=size(img);
fgImg=zeros(sz,'uint8');    %08/19/2021, changed from uint16 to unit8
ndim=ndims(img);

for i=1:size(pts,1)
    bndBox=getPointBndBox(pts(i,:), cutwidth, 0, sz, 'small');
    ymin=bndBox(1,1);
    ymax=bndBox(2,1);
    xmin=bndBox(1,2);
    xmax=bndBox(2,2);
%     ymin=max([pts(i,1)-cutwidth_xy, 2]);
%     xmin=max([pts(i,2)-cutwidth_xy,2]);
%     ymax=min([pts(i,1)+cutwidth_xy, nY-1]);
%     xmax=min([pts(i,2)+cutwidth_xy, nX-1]);
    if ndim==2
        fgImg(ymin:ymax, xmin:xmax)=fgImg(ymin:ymax, xmin:xmax)+1;
    else
        zmin=bndBox(1,3);
        zmax=bndBox(2,3);
%         zmin=pts(i,3)-cutwidth_z; zmin=max(1,min(zmin,nZ));
%         zmax=pts(i,3)+cutwidth_z; zmax=max(1, min(zmax,nZ));
        fgImg(ymin:ymax, xmin:xmax, zmin:zmax)=fgImg(ymin:ymax, xmin:xmax, zmin:zmax)+1;
    end
end
end