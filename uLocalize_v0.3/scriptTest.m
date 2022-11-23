img=imread('C:\Users\Bin\OneDrive - Johns Hopkins University\Programs\Matlab\TransTrack\TransTrackTest\ATF4_3.tif');
p=initAirLocalizePara();
p.sigma_xy=0.8008;
smooth=bpass_filter_2D_fourier(img,p.filter.nlo,p.sigma_xy,p.filter.nhi,p.filter.width); 
smooth=smooth-mean(smooth(:));
% msk=drawMask(img);
pts=predetectWithThreshInSTD(smooth, 6, 1, msk);
bgImg=genBGImageByRemovePts(img,pts,3);
cutwidth=3;
thickness=1;
ROISizeOption='small';
final_pix=zeros(size(pts,1),4);
for i=1:size(pts,1)
    [imgCorr, pt2, bndBox, bgCrop]=genBGCorrByLinInterp2D(img, bgImg, pts(i,:), cutwidth, thickness, ROISizeOption);
    %The original airlocalize program
%     [y0,x0,N0,~,~,chi2]=gaussian_2Dmask_small2_chi2(imgCorr, pt2(1), pt2(2), p.cutsize, p.sigma_xy,1,1,p.maxcount,p.tol);
%     final_pix(i,1)=y0(end)+bndBox(1,1)-0.5;
%     final_pix(i,2)=x0(end)+bndBox(1,2)-0.5;
%     final_pix(i,3)=N0(end);
%     final_pix(i,4)=chi2;

%     %The Gaussian Mask Algorithm
%     pt=gaussianMask2D(imgCorr, pt2, p.sigma_xy, p.maxcount, p.tol);
%     final_pix(i,1)=pt(1)+bndBox(1,1)-0.5;
%     final_pix(i,2)=pt(2)+bndBox(1,2)-0.5;
%     final_pix(i,3)=pt(3);
%     final_pix(i,4)=pt(4);
    
    %The Gaussian Fit algorithm
    pt=gaussianFit2D(imgCorr,pt2,'integrated gaussian');
%     pt=gaussianFit2D(imgCorr,pt2,'standard gaussian');
    final_pix(i,1)=pt(1)+bndBox(1,1)-0.5;
    final_pix(i,2)=pt(2)+bndBox(1,2)-0.5;
    final_pix(i,3)=pt(3);
    final_pix(i,4)=pt(4);    
end
imshow(img,[]);
hold on;
plot(final_pix(:,2),final_pix(:,1),'sq');

