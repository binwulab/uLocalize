function data=getBGDataAroundPoint(img,fgImg,pt,cutwidth,thickness)
%Objective: generate a region in the image around a point. Will use it to
% calculate bg correction for the objecti
%Input
%   img: the original image, 2D or 3D
%   fgImg: the foreground image with all candiates points and their immediate surrounding marked >0 
%   pt: the candidate point: [Y,X,Z,Int], or [Y,X,Int] for 2D
%   cutwidth: the size of the image: [cutwidth_xy, cutwidth_z]
%           image returned: [2*cutwidth+1, 2*cutwidth+1]
%   thickness: the thickness of around the region to be cut, used for bg correction
%Output
%   data: the background data [y,x,z,int] or [y,x,int] used to fit a plane
%History
%   BW, Feb 2021

sz=size(img);
numdim=ndims(img);
%define the boxes: [y1,x1; y4,x4] and [y2,x2; y3,x3]
bndBox14=getPointBndBox(pt, cutwidth,thickness,sz,'large');
y1=bndBox14(1,1);
x1=bndBox14(1,2);
y4=bndBox14(2,1);
x4=bndBox14(2,2);

bndBox23=getPointBndBox(pt, cutwidth,thickness,sz,'small');
y2=bndBox23(1,1);
x2=bndBox23(1,2);
y3=bndBox23(2,1);
x3=bndBox23(2,2);

% if y1==1; y2=2; end     %If it is at the boundary, keep only one line
% if x1==1; x2=2; end
% if y4==ny; y3=ny-1; end
% if x4==nx; x3=nx-1; end

%Now change y2,x2,y3 and x3 to relative coordinate
y2=max([y2-y1+1,1]); 
x2=max([x2-x1+1,1]);
y3=min([y3-y1+1,y4-y1+1]);
x3=min([x3-x1+1,x4-x1+1]);

if numdim == 3
    z1=bndBox14(1,3);
    z4=bndBox14(2,3);
    
    z2=bndBox23(1,3);
    z3=bndBox23(2,3);
    %change z2, z3, to relative coordinates
    z2=max([z2-z1+1,1]);
    z3=min([z3-z1+1,z4-z1+1]);
    [X,Y,Z]=meshgrid(x1:x4, y1:y4, z1:z4);
    X(y2:y3, x2:x3, z2:z3)=nan;    %make the center of the box NAN
    fgImg_crop=fgImg(y1:y4,x1:x4,z1:z4)-1; %The foreground is uint8 containing number of obj in each pixel, this will remove the current points and its surroundings. if there is overlapping points, it will not be used
    img_crop=img(y1:y4,x1:x4,z1:z4);
    indx=find(~isnan(X) & (fgImg_crop < 1));   %Choosing only the boundary and zero foreground pixels (bg). Note somehow the nan pixel is set to 0.
    data=zeros(numel(indx),4);  %The data used for fitting the plane
    data(:,1)=Y(indx)-0.5;
    data(:,2)=X(indx)-0.5;
    data(:,3)=Z(indx)-0.5;
    data(:,4)=img_crop(indx);
else
    [X,Y]=meshgrid(x1:x4, y1:y4);
    X(y2:y3, x2:x3)=nan;    %make the center of the box NAN
    fgImg_crop=fgImg(y1:y4,x1:x4)-1; %The foreground is uint16 containing number of obj
    img_crop=img(y1:y4,x1:x4);
    indx=find(~isnan(X) & (fgImg_crop < 1));   %Choosing only the boundary  and zero foreground pixels (bg). Note somehow the nan pixel is set to 0.
    data=zeros(numel(indx),3);  %The data used for fitting the plane
    data(:,1)=Y(indx)-0.5;  %Note: the center of a pixel
    data(:,2)=X(indx)-0.5;
    data(:,3)=img_crop(indx); 
end
%now keeping 80% of the pixels (removing outliers)
nd=numel(indx);
if nd>10
    [~,indx]=sort(data(:,4),'ascend');
    data=data(indx,:);
    n1=round(0.1*nd);
    n2=round(0.9*nd);
    data=data(n1:n2,:);
end
end