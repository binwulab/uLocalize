function [pt,status]=gaussianMaskWrapper(Im, pt0, sigma, maxIter, tol, bw)
numdim=ndims(Im);
if numdim==3
    [pt,status]=gaussianMask3D(Im, pt0, sigma, maxIter, tol, bw);
else
    [pt,status]=gaussianMask2D(Im, pt0, sigma(1), maxIter, tol, bw);
end
end