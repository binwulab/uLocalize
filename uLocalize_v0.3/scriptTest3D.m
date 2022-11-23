img=tiffread5('C:\Users\Bin\OneDrive - Johns Hopkins University\Programs\Matlab\Airlocalize\AIRLOCALIZE1_stable_2012_07_26\AIRLOCALIZE1_stable\bwExtension\C2-A60_NT_1.tif',1,31);
p=initAirLocalizePara();
p.sigma_xy=1.1918;
p.sigma_z=0.9246;
smooth=double(bpass_filter_3D_fast_memory_safe3(img,p.filter.nlo,p.sigma_xy,p.sigma_z,p.filter.nhi,p.filter.width,p.filter.numdim)); 
smooth=smooth-mean(smooth(:));
% msk=drawMask(img);
cutwidth=[3,3];
thickness=1;
ROISizeOption='small';
pts=predetectWithThreshInSTD(smooth, 8, 1, msk);
bgImg=genBGImageByRemovePts(img,pts,cutwidth);
final_pix=zeros(size(pts,1),7);
for i=1:size(pts,1)
    [imgCorr, pt2, bndBox, bgCrop]=genBGCorrByLinInterp3D(img, bgImg, pts(i,:), cutwidth, thickness, ROISizeOption);
%     [y0,x0,z0, N0,err0]=gaussian_3Dmask_small4(imgCorr, pt2(1), pt2(2), pt2(3), p.cutsize, p.sigma_xy,p.sigma_z,1,1,1,p.maxcount,p.tol);
%     final_pix(i,1)=y0(end)+bndBox(1,1);
%     final_pix(i,2)=x0(end)+bndBox(1,2);
%     final_pix(i,3)=z0(end)+bndBox(1,3);
%     final_pix(i,4)=N0(end);
%     final_pix(i,5)=err0;

%     %The Gaussian Mask Algorithm
    pt=gaussianMask3D(imgCorr, pt2, [p.sigma_xy,p.sigma_z],p.maxcount,p.tol);
    final_pix(i,1)=pt(1)+bndBox(1,1);
    final_pix(i,2)=pt(2)+bndBox(1,2);
    final_pix(i,3)=pt(3)+bndBox(1,3);
    final_pix(i,4)=pt(4);   %N0
    final_pix(i,5)=pt(5);   %err

    %The Gaussian Fit algorithm
%     pt=gaussianFit3D(imgCorr, pt2, 'integrated gaussian');
%     final_pix(i,1)=pt(1)+bndBox(1,1)-0.5;
%     final_pix(i,2)=pt(2)+bndBox(1,2)-0.5;
%     final_pix(i,3)=pt(3)+bndBox(1,3)-0.5;
%     final_pix(i,4)=pt(4); %N0
%     final_pix(i,5)=pt(5); %sigma_xy
%     final_pix(i,6)=pt(6); %sigma_z
end
imshow(max(img,[],3),[600,5000]);
hold on;
plot(final_pix(:,2)-0.5,final_pix(:,1)-0.5,'sq');

