function [data,bndBox,pixelIdxList2,bw2]=getBGDataAroundObj(img,bw,pixelIdxList,thickness,bg_ext)
%Objective: generate a region in the image around an object. Will use it to
% calculate bg correction for the objecti
%Input:
%   img: the original image, 2D or 3D
%   bw: a mask file containing all objects in the image, 0 would be bg
%   pixelIdxList: the current object, represented as a linear indices in img
%   thickness: the thickness of the background region
%   bg_ext: the number of pixel extended from the boundary of the object
%Output
%   data: the background data [y,x,z,int] or [y,x,int] used to fit a plane
%   bndBox: the bounding box for the cropped image [ymin,xmin,zmin; ymax,xmax,zmax]
%   pixelIdxList: the linear indices in the cropped image
%   bw2: the cropped mask image. Note it may contain other obj
%History
%   BW, Jan 2021

numdim=ndims(img);
imSize=size(img);
bg_ext=round(bg_ext);
if thickness > bg_ext   %So thickness cannot exceed bg_ext. 
    thickness = bg_ext;
end
[bndBox, pixelIdxList2, bw2]=getObjBndBox(pixelIdxList, imSize, bg_ext, bw);
%The outer box
y1=bndBox(1,1);
y4=bndBox(2,1);
x1=bndBox(1,2);
x4=bndBox(2,2);
%The inner box
y2=y1+thickness;
x2=x1+thickness;
y3=y4-thickness;
x3=x4-thickness;

%change to relative coordinates
y2=max([y2-y1+1,1]);
x2=max([x2-x1+1,1]);
y3=min([y3-y1+1,y4-y1+1]);
x3=min([x3-x1+1,x4-x1+1]);

if numdim == 3
    z1=bndBox(1,3);
    z4=bndBox(2,3);
    
    z2=z1+thickness;
    z3=z4-thickness;
    
    %change to relative coordinates
    z2=max([z2-z1+1,1]);
    z3=min([z3-z1+1,z4-z1+1]);
    
    [X,Y,Z]=meshgrid(x1:x4,y1:y4,z1:z4);
    X(y2:y3,x2:x3,z2:z3)=nan;   %set the middle of the cropped coord to nan
    indx=find(~isnan(X) & bw2==0);  %take out the surrounding background pixels 
    data=zeros(numel(indx),4);  %The data used for fitting the plane
    data(:,1)=Y(indx)-0.5;
    data(:,2)=X(indx)-0.5;
    data(:,3)=Z(indx)-0.5;
    data(:,4)=img(indx);
    %now keeping only half of the pixels (removing outliers)
    [~,indx]=sort(data(:,4),'ascend');  
    data=data(indx,:);  %sort according to the intensity
elseif numdim==2
    [X,Y]=meshgrid(x1:x4,y1:y4);
    X(y2:y3,x2:x3)=nan;
    indx=find(~isnan(X) & bw2==0);
    data=zeros(numel(indx),3);  %The data used for fitting the plane
    data(:,1)=Y(indx)-0.5;
    data(:,2)=X(indx)-0.5;
    data(:,3)=img(indx);
    %now keeping only half of the pixels (removing outliers)
    [~,indx]=sort(data(:,3),'ascend');
    data=data(indx,:);  %sort according to the intensity
end

%keep only the middle range of the data
nd=numel(indx);
if nd>4
    n1=round(nd/6);
    n2=round(5*nd/6);
    data=data(n1:n2,:);
end
end