function plotGranule(cc)
imSize=cc.ImageSize;
numdim=numel(imSize);
for i=1:cc.NumObjects
    if numdim==3
        [yobj, xobj, zobj]=ind2sub(imSize, cc.PixelIdxList{i});
        plot(xobj, yobj,'y.');
    end
end
end