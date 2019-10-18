% function MibiCompareColorinTwoImages

impath1 = 'Point1/TifsNoNoise/dsdna (r) hla1 (g) (RGB).tif';
impath2 = 'Point1/TifsNoNoise/dsDNA (r) ki67 (g) (RGB).tif';

sumColorVec1 = MibiQuantifyColorInImage(impath1);
sumColorVec2 = MibiQuantifyColorInImage(impath2);

% get ratio of yellow and green
ratio1 = sumColorVec1(4)/(sumColorVec1(4)+sumColorVec1(2));
ratio2 = sumColorVec2(4)/(sumColorVec2(4)+sumColorVec2(2));

% plot
figure;
bar([ratio1,ratio2],'FaceColor',[0.5,0.5,0.5]);
set(gca,'xTickLabel',{'HLA-1','Ki67'});