% check new and old segmentation
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path1 = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
path2 = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
pointNumber=2;

load([path1,'/Point',num2str(pointNumber),'/dataDeNoiseCohortNoAgg.mat']);
old=load([path1,'/Point',num2str(pointNumber),'/segmentationParams.mat']);
new=load([path2,'/Point',num2str(pointNumber),'/segmentationParams.mat']);

figure; imagesc(old.cellPerimNewMod);
plotbrowser on;
figure; imagesc(new.cellPerimNewMod);
plotbrowser on;

maxv=20;
rgb_image = MibiGetRGBimageFromMat(countsNoNoiseNoAgg(:,:,7),maxv);
rgb_image_perim1= imoverlay(rgb_image , old.cellPerimNewMod , [1 .3 .3]);
figure; imagesc(rgb_image_perim1);
plotbrowser on;

rgb_image_perim2= imoverlay(rgb_image , new.cellPerimNewMod , [1 .3 .3]);
figure; imagesc(rgb_image_perim2);
plotbrowser on;

