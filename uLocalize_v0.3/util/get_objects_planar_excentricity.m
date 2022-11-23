function [stats,cc2] = get_objects_planar_excentricity(cc)

if cc.NumObjects == 0
    stats = 0;
    cc2 = cc;
    return
end
% with a 3D stack, computing the 2D projections of 3D objects 
if size(cc.ImageSize,2) == 3
    cc2 = cc;
    cc2.ImageSize = cc.ImageSize(1:2);
    cc2.Connectivity = 8;
    cc2 = rmfield(cc2,'PixelIdxList');
    for i=1:length(cc.PixelIdxList)

        %retrieving the (x,y) positions
        PixList = mod( cc.PixelIdxList{1,i}, cc.ImageSize(1)*cc.ImageSize(2) ); %z projection (using linear indexing)
        PixList = unique(PixList); %removing non unique positions
        PixList(PixList == 0) = cc.ImageSize(1)*cc.ImageSize(2) ;
        [x,y,~] = ind2sub(cc.ImageSize(1:2), PixList ); 
        clear('PixList');

        %replacing the projected region in the 2D set of objects
        cc2.PixelIdxList{1,i} = sub2ind(cc2.ImageSize,x,y);
        clear('x','y');
    end
    stats = regionprops(cc2,'Eccentricity');
else
    stats = regionprops(cc,'Eccentricity');
end

stats = [stats.Eccentricity];


end