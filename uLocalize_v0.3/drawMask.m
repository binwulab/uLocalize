function mask=drawMask(img)
sz=size(img);
if(numel(sz)==3)    %max projection
    img=max(img,[],3);
end
h=imshow(img,[1.1*min(img(:)), 0.9*max(img(:))]);
hold on;
hroi=impoly();
mask=createMask(hroi,h);
end