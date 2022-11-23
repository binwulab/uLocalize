%% 3D Fish image
% img=tiffread5('C:\Users\Bin\OneDrive - Johns Hopkins University\data\Test\FISH\MEF\C1-steady_crop.tif',1,21);
% % img=tiffread5('E:\OD\OneDrive - Johns Hopkins University\data\Test\FISH\MEF\C1-steady_crop.tif',1,21);
% img=tiffread5('C:\Users\bwu20\Desktop\Lauren PBody\C1-mbsmef_16_steadystate_dox_xy01.tif',1,21);
img=tiffread5('C:\Users\bwu20\Desktop\Lauren PBody\C4-steady.tif',1,21);    %P-body channel
%%
%% detection with connected component: dense mode
p=uLocalizeInitPara('numdim',3, 'sigma_xy', 1.38, 'sigma_z',1.05, 'numDilation',3, 'aMin', 20, 'eMax', 0.9);
maskType=1;
% mask=drawMask(img);
mask=ones(size(img(:,:,1)));
p.detectionMode='granule';
p.filterMethod='LOGRAJ';
p.filter.filterSigma=p.sigma_xy*2;
% smooth=uLocalizeFilter(img,p);
%%
thresh=10;
p.bg_extension=3;
p.thickness=1;
%cc=predetectObjectsWithThreshInSTD(smooth, thresh, maskType,mask,[p.aMin, p.aMax],[p.eMin,p.eMax], true);
[final_pts, p, ~, detection_res]=uLocalizeDetection(img,p,maskType,mask,'thresh',thresh,'method','GaussianMask','calibration',true);
imshow(max(img,[],3),[400,15000]);
hold on;
plot(final_pts(:,2), final_pts(:,1),'ro')
plotGranule(detection_res);
%%
[final_pts, p, detection_res]=uLocalizeDetection(img,p,maskType,mask,'thresh',thresh,'method','GaussianMask','calibration',true);
spotSingle=detection_res.single;
spotMultiple=detection_res.multiple;
imshow(max(img,[],3),[400,1500]);
hold on;
plot(spotSingle.fitRes(:,2),spotSingle.fitRes(:,1),'gs');
% plot(spotMultiple.fitRes(:,2),spotMultiple.fitRes(:,1),'bo');
cc=detection_res.cc;
plot(spotMultiple.fitRes(1:cc.nObjects(2),2),spotMultiple.fitRes(1:cc.nObjects(2),1),'bo');
indx=find(spotMultiple.nRNA>=1);
plot(spotMultiple.fitRes(indx,2),spotMultiple.fitRes(indx,1),'ro');
sum(spotMultiple.nRNA)

%% detection with local maximum
p=uLocalizeInitPara('numdim',3, 'sigma_xy', 1.38, 'sigma_z',1.05, 'cutsize',2);
maskType=1;
% mask=drawMask(img);
% mask=[];
% mask=drawMask(img);
thresh=15;
[final_pt]=uLocalizeDetection(img,p,maskType,mask,'thresh',thresh,'method','GaussianMask');
imshow(max(img,[],3),[400,1000]);
hold on;
plot(final_pt(:,2),final_pt(:,1),'gs');

%% 2D image
img=imread('C:\Users\Bin\OneDrive - Johns Hopkins University\Programs\Matlab\uLocalize\ATF4_3.tif');
% detection
maskType=0;
mask=[];
% mask=drawMask(img);
p=uLocalizeInitPara('numdim',2,'sigma_xy', 0.8008, 'numDilation',3, 'aMin', 3, 'eMax', 0.7);
thresh=4;
[spotSingle, spotMultiple, p, cc]=uLocalizeDetectionCC(img,p,maskType,mask,'thresh',thresh,'method','GaussianMask','calibration',true);
imshow(img,[]);
plot(spotSingle.fitRes(:,2),spotSingle.fitRes(:,1),'gs');
% plot(spotMultiple.fitRes(:,2),spotMultiple.fitRes(:,1),'bo');
plot(spotMultiple.fitRes(1:cc.nObjects(2),2),spotMultiple.fitRes(1:cc.nObjects(2),1),'bo');
indx=find(spotMultiple.nRNA>=1);
plot(spotMultiple.fitRes(indx,2),spotMultiple.fitRes(indx,1),'ro');
sum(spotMultiple.nRNA)
