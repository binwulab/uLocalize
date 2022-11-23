function bndBox=getPointBndBox(pt,cutwidth,thickness,imsize,ROISizeOption)
%Objective: Given a point, calculate the bounding box for cropping image
%Input
%   pt: the point of interest, [Y, X, Z]
%   cutwidth: the size of ROI in pixel [cutwidth_xy, cutwidth_z]
%   thickness: the surrounding ring thickness to calculate the background
%   imsize: the size of the image [nY, nX, nZ], determining dimension of
%           results
%   ROISizeOption: 'large': including the surrounding ring with thickness
%           'small', only including the region with cutwidth
%Output
%   bndBox: [ymin, xmin, zmin; ymax, xmax, zmax] for 3D
%           [ymin, xmin; ymax zmax] for 2D 
yc=ceil(pt(1));
xc=ceil(pt(2));
cutwidth_xy=floor(cutwidth(1));
ny=imsize(1);
nx=imsize(2);
numdim=numel(imsize);       %Note: the problem dim is determined by imsize
if numdim == 3
    zc=ceil(pt(3));
    cutwidth_z=floor(cutwidth(2));
    nz=imsize(3);
end
if strcmpi(ROISizeOption,'large')
    ymin=max([yc-cutwidth_xy-thickness,1]);
    xmin=max([xc-cutwidth_xy-thickness,1]);
    ymax=min([yc+cutwidth_xy+thickness,ny]);
    xmax=min([xc+cutwidth_xy+thickness,nx]);
    if numdim ==3
        zmin=max([zc-cutwidth_z-thickness,1]);
        zmax=min([zc+cutwidth_z+thickness,nz]);
    end
else
    ymin=max([yc-cutwidth_xy,1]);
    xmin=max([xc-cutwidth_xy,1]);
    ymax=min([yc+cutwidth_xy,ny]);
    xmax=min([xc+cutwidth_xy,nx]);
    if numdim ==3
        zmin=max([zc-cutwidth_z,1]);
        zmax=min([zc+cutwidth_z,nz]);
    end
end
if numdim==3
    bndBox=[ymin,xmin,zmin; ymax,xmax,zmax];
else
    bndBox=[ymin,xmin; ymax,xmax];
end
end