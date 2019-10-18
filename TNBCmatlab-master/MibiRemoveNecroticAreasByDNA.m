% Mibi get necrosis mask from dna channel and remove it from dna to assist
% segmentation


path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/';
pointNumber = 15;
channelInd=7;

load([path,'/Point',num2str(pointNumber),'/dataDeNoise.mat']);
data = countsNoBg(:,:,channelInd);

%%
cap = 200;
t1=0.008;
t2=0.038;
removeVal = 10;
gausRad1=1;
gausRad2=1;
tsize=10000;

MibiGausCapAndPlot(data,25,1); plotbrowser on;

mask1 = MibiGetMask(data,cap,t1,gausRad1);
mask2 = MibiGetMask(data,cap,t2,gausRad2);
necrosisMask = mask1-mask2;
figure; imagesc(necrosisMask); plotbrowser on;

se = strel('disk',5);
necrosisMask1 = imclose(necrosisMask,se);
figure; imagesc(necrosisMask1); plotbrowser on;
necrosisMask2 = imopen(necrosisMask1,se);
figure; imagesc(necrosisMask2); plotbrowser on;
% remove elements with size < t
necrosisMask3 = necrosisMask2;
CC = bwconncomp(necrosisMask2);
stats = regionprops(CC,'Area');
for i=1:length(stats)
    if (stats(i).Area < tsize)
        coord = CC.PixelIdxList{i};
        necrosisMask3(coord) = 0;
    end
end
figure; imagesc(necrosisMask3); plotbrowser on;
dataNew = MibiRemoveBackgroundByMaskSingleChannel(countsNoBg(:,:,channelInd),logical(necrosisMask3),removeVal);
MibiGausCapAndPlot(dataNew,25,1); plotbrowser on;
countsNoBg(:,:,channelInd) = dataNew;
save([path,'/Point',num2str(pointNumber),'/dataDeNoise.mat'],'totalIonFiltSum','countsAllSFiltCRSum','countsNoNoise','massDS','IntNormD','pointNumber','countsNoBg');