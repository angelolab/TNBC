% Mibi compare Halo segmentation to mine

PointNumber = 3;
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
load([path,'/Point',num2str(PointNumber),'/dataDeNoiseCohortNoAgg.mat']);
load([path,'/Point',num2str(PointNumber),'/segmentationParams.mat']);
%haloDataRGB = imread([path,'/HALO/ExportFromHALO/Halo archive 2017-07-31 09-59/Images/NucleiGaus_p3_6_job18_MarkupActual.tif']);
%haloDataRGB = imread([path,'/HALO/ExportFromHALO/Halo archive 2017-07-31 11-45/Images/Nuclei_p3_5_job13_MarkupActual.tif']);
haloDataRGB = imread([path,'/HALO/ExportFromHALO/Halo archive 2017-07-31 13-21/Images/Nuclei_p3_5_job19_MarkupActual.tif']);

haloData = rgb2gray(haloDataRGB);
figure; imagesc(haloData);

% overlay my mask on halo segmentation
rgb_image_perim= imoverlay(haloData , cellPerimNewMod , [1 .3 .3]);
figure; imagesc(rgb_image_perim);

% plot nuclei
parula();
rgb_image_nuc = MibiGetRGBimageFromMat(countsNoNoiseNoAgg(:,:,7),15);
rgb_image_perim_nuc= imoverlay(rgb_image_nuc , cellPerimNewMod , [1 .3 .3]);
figure; imagesc(rgb_image_perim_nuc);
