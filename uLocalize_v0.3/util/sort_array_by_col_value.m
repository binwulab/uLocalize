function [sort_arr,indices] = sort_array_by_col_value(arr,ncol,orientation)
%arr is the array, 
%col is the column which values decide of the ordering of
%the rows
%orientation is a string, either 'ascend' or 'descend'

ny = size(arr,2);
if ny<ncol
    disp('your array does not have enough columns - get yourslef a greek temple');
end

[sortVal indices] = sort(arr(:,ncol),orientation);
nx = size(arr,1);
sort_arr(1:nx,1:ny) = arr(indices(1:nx),1:ny);
clear('nx','ny','sortVal');
end