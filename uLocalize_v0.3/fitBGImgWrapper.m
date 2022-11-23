function [para, bg_sd, res_mean]=fitBGImgWrapper(data, bgCorrMethod)
%Obj: given the data in the format [x,y,int] for 2D and [x,y,z,int] for 3D,
%   fit the data to a function (hyperplane). It is used to calculate the bg image
%   surrounding an object, 
%Input: 
%   data: [x,y,int] for 2D, [x,y,z,int] for 3D
%   fitMethod: method to be used for the background 
%   fitPara: the parameter passed for a given method, not used for now
%Output:
%   para: the parameter for the linear fit to the data
%       To use the para to generate a 4D plane:
%       para(1)*(Y-0.5)+para(2)*(X-0.5)+para(3)*(Z-0.5)+para(4) 
%       3D plane
%       para(1)*(Y-0.5)+para(2)*(X-0.5)+para(3)
%   bg_sd: the std of the difference betwwen the fitted value and the data
%   res_mean: the mean deviation between the fitted value and the data
%History
%   BW, Aug 2021
if nargin<2
    bgCorrMethod = 'plane';
end
numdim=size(data,2)-1;  %Assume the dimension of the problem is given in the data format
if numdim == 3
    switch lower(bgCorrMethod)
        case {'plane', 'local plane'}
            para=fitBGImgW4DPlane(data);
            Ires=data(:,4)-(para(1)*data(:,1)+para(2)*data(:,2)+para(3)*data(:,3)+para(4));
        case {'quad','local quad'}
            para=fitBGImgW4DQuad(data);
            Ires=data(:,4)-(para(1)*data(:,1)+para(2)*data(:,2)+para(3)*data(:,3)+para(4)+para(5)*data(:,3).^2);
        case {'median','local median'}
            para = median(data(:,4));
            Ires=data(:,4)-para;
        otherwise
            disp('fitBGImgWrapper 3D: no such background correction method, dear');
            para=[];
            bg_sd=[];
            res_mean=[];
            return;
    end
else
    switch lower(bgCorrMethod)
        case {'plane', 'local plane'}
            para=fitBGImgW3DPlane(data);
                Ires=data(:,3)-(para(1)*data(:,1)+para(2)*data(:,2)+para(3));
        case {'median','local median'}
            para = median(data(:,3));
            Ires=data(:,3)-para;
        otherwise
            disp('fitBGImgWrapper 2D: no such background correction method, dear');
            para=[];
            bg_sd=[];
            res_mean=[];
            return;  
    end
end
if nargout>1
    bg_sd=std(Ires);
end
if nargout > 2
    res_mean=mean(Ires);
end
end



function par=fitBGImgW3DPlane(data)
%Fit the intensity I vs (y,x) to a 3D plane with equation I=a*y+b*x+c
%Input: 
%   data: (y,x,int)
%Output
%   para: [a,b,c]
data=double(data);
npts=size(data,1);
s1 = sum(data(:,1));
s2 = sum(data(:,2));
s3 = sum(data(:,3));
s11 = sum( data(:,1).*data(:,1) );
s22 = sum( data(:,2).*data(:,2) );
s12 = sum( data(:,1).*data(:,2) );
s31 = sum( data(:,3).*data(:,1) );
s32 = sum( data(:,3).*data(:,2) );

fitmat = [s11 s12 s1; s12 s22 s2; s1 s2 npts];
DD = det(fitmat);

%if matrix inversion is impossible, I just return the average intensity as
%the result
if (abs(DD)<1e-6)
    disp('Warning fitBGImgW3Dplane: the matrix determinant is 0');
    par = [0;0;mean(data(:,3))];
else
    v = [s31;s32;s3];
    par = fitmat\v;
end
end

function par=fitBGImgW4DPlane(data)
%Fit the intensity I vs (y,x) to a 3D plane with equation I=a*y+b*x+c*z+d
%The solution is a linear least square problem
%Input: 
%   data: (y,x,z,int)
%Output
%   para: [a,b,c,d]
data=double(data);
npts=size(data,1);
s1 = sum(data(:,1));
s2 = sum(data(:,2));
s3 = sum(data(:,3));
s4 = sum(data(:,4));
s11 = sum( data(:,1).*data(:,1) );
s22 = sum( data(:,2).*data(:,2) );
s33 = sum( data(:,3).*data(:,3) );
s12 = sum( data(:,1).*data(:,2) );
s13 = sum( data(:,1).*data(:,3) );
s23 = sum( data(:,2).*data(:,3) );
s41 = sum( data(:,4).*data(:,1) );
s42 = sum( data(:,4).*data(:,2) );
s43 = sum( data(:,4).*data(:,3) );
fitmat = [s11 s12 s13 s1; s12 s22 s23 s2; s13 s23 s33 s3; s1 s2 s3 npts];
DD = det(fitmat);

%if matrix inversion is impossible, I just return the average intensity as
%the result
if (abs(DD)<1e-6)
    disp('Warning fitBGImgW4DPlane: the matrix determinant is 0');
    par = [0;0;0;mean(data(:,4))];
else
    v = [s41;s42;s43;s4];
    par = fitmat\v;
end
end

function par=fitBGImgW4DQuad(data)
%Fit the intensity I vs (y,x) to a 3D plane with equation I=a*y+b*x+c*z+d+e*z^2
%The solution is a linear least square problem
%Input: 
%   data: (y,x,z,int)
%Output
%   para: [a,b,c,d,e]
data=double(data);
npts=size(data,1);
s1 = sum(data(:,1));    %Y
s2 = sum(data(:,2));    %X
s3 = sum(data(:,3));    %Z
s4 = sum(data(:,4));    %I
s11 = sum( data(:,1).*data(:,1) );  %YY
s22 = sum( data(:,2).*data(:,2) );  %XX
s33 = sum( data(:,3).*data(:,3) );  %ZZ
s12 = sum( data(:,1).*data(:,2) );  %YX
s13 = sum( data(:,1).*data(:,3) );  %YZ
s23 = sum( data(:,2).*data(:,3) );  %XZ
s41 = sum( data(:,4).*data(:,1) );  %IY
s42 = sum( data(:,4).*data(:,2) );  %IX
s43 = sum( data(:,4).*data(:,3) );  %IZ
s1_ZZ = sum(data(:,1).*data(:,3).^2);   %YZ^2
s2_ZZ =  sum(data(:,2).*data(:,3).^2);  %XZ^2
s3_ZZ =  sum(data(:,3).^3);             %Z^3
s4_ZZ =  sum(data(:,4).*data(:,3).^2);   %IZ^2


fitmat =   [s11     s12     s13     s1  s1_ZZ; ...
            s12     s22     s23     s2  s2_ZZ; ...
            s13     s23     s33     s3  s3_ZZ; ...
            s1      s2      s3      npts s33; ...
            s1_ZZ   s2_ZZ   s3_ZZ   s33  s3_ZZ];
DD = det(fitmat);

%if matrix inversion is impossible, I just return the average intensity as
%the result
if (abs(DD)<1e-6)
    disp('Warning fitBGImgW4DQuad: the matrix determinant is 0');
    par = [0;0;0;mean(data(:,4));0];
else
    v = [s41;s42;s43;s4;s4_ZZ];
    par = fitmat\v;
end
end

