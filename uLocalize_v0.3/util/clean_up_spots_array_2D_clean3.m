function [fspots2,fpix2,n_wrong,n_double] = ...
    clean_up_spots_array_2D_clean3(fpix,vox_size,im_size,n_exclusion,numdim)
%series of checks on the detected spots to remove aberrant detections
if numdim == 3
    nx = im_size(1);
    ny = im_size(2);
    nz = im_size(3);
    Icolnum = 4;
else
    nx = im_size(1);
    ny = im_size(2);
    Icolnum = 3;
end

fpix = sort_array_by_col_value(fpix,Icolnum,'descend');
%% removing the rows corresponding to iterations of the gaussian mask that did not converge 
%(signaled in the gaussian mask algorithm by [x y z N] = [-1 -1 -1 -1])

not_wrong = logical((fpix(:,Icolnum) ~= -1).*(fpix(:,1) >= 0).*(fpix(:,1) < nx).*(fpix(:,2) >= 0 ).*(fpix(:,2) < ny ));
if numdim == 3
    not_wrong = not_wrong.*( fpix(:,3) >= 0 ).*( fpix(:,3) < nz);
end

n_wrong = size(fpix,1) - sum(not_wrong);
% save('D:\junk\cleanup.mat','not_wrong');
fpix = fpix(not_wrong,:);
clear('not_wrong');

%% removing duplicate spots that converged within n_exclusion pixel from each other

%I create a stack in which each spot is sent to its pixel, with its value
%corresponding to its intensity rank (1: least intense spot)
tmp = [uint16(fpix(:,1:numdim)),uint16((size(fpix,1):-1:1)')];
if numdim ==3
    spots_stack = convert_coordinates_to_stack2(tmp,nx,ny,nz);
else
    spots_stack = convert_coordinates_to_stack2(tmp,nx,ny);
end

%I check whether each spot is the maximal ranking within an ROI of
%size n_exclusion, otherwise I discard it
maxima = remove_non_localmax_pixels(spots_stack,tmp,n_exclusion);

n_double = size(fpix,1) - size(maxima,1);
fpix2 = fpix(size(fpix,1)-maxima(:,4)+1,:);

fspots2 = fpix2;
fspots2(:,1) = fpix2(:,1)*vox_size(1);
fspots2(:,2) = fpix2(:,2)*vox_size(2);
if numdim == 3
    fspots2(:,3) = fpix2(:,3)*vox_size(3);
end

%% the end
clear('spots_stack','tmp','minima','fpix','not_wrong','t','indices');

end

function maxima = remove_non_localmax_pixels(stack,pts,ROIsize)
    %removes a point if there is another one with a higher intensity
    %within 'ROIsize' pixels in any direction (faster)

if ndims(stack) == 3
    [nx,ny,nz] = size(stack);
    maxima = zeros(size(pts,1),4);
    
    j=1;
    for i = 1:size(pts,1)
        
        %generate ROI around the current point
        xmin = pts(i,1)-ROIsize; xmin = max(1,min(xmin,nx));
        xmax = pts(i,1)+ROIsize; xmax = max(1,min(xmax,nx));
        ymin = pts(i,2)-ROIsize; ymin = max(1,min(ymin,ny));
        ymax = pts(i,2)+ROIsize; ymax = max(1,min(ymax,ny));
        zmin = pts(i,3)-ROIsize; zmin = max(1,min(zmin,nz));
        zmax = pts(i,3)+ROIsize; zmax = max(1,min(zmax,nz));
        ROI = stack(xmin:xmax,ymin:ymax,zmin:zmax);
        
        if max(max(max(ROI))) == pts(i,4)
            maxima(j,:) = pts(i,:);
            j=j+1;
        end
    end
    if j>1
        maxima = maxima(1:j-1,:);     
    end   
elseif ndims(stack)==2
    [nx,ny] = size(stack);
    maxima = zeros(size(pts,1),3);
    
    j=1;
    for i = 1:size(pts,1)
        
        %generate ROI around the current point
        xmin = pts(i,1)-ROIsize; xmin = max(1,min(xmin,nx));
        xmax = pts(i,1)+ROIsize; xmax = max(1,min(xmax,nx));
        ymin = pts(i,2)-ROIsize; ymin = max(1,min(ymin,ny));
        ymax = pts(i,2)+ROIsize; ymax = max(1,min(ymax,ny));
        ROI = stack(xmin:xmax,ymin:ymax);
        
        if max(max(max(ROI))) == pts(i,3)
            maxima(j,:) = pts(i,:);
            j=j+1;
        end
    end
    if j>1
        maxima = maxima(1:j-1,:);     
    end
end


clear('nx','ny','nz','i','j','xmin','xmax','ymin','ymax','zmin','zmax','ROI');
end