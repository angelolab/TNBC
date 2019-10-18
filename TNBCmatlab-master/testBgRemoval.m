mask = MibiGetMask(countsAllSFiltCRSum(:,:,10),10);
countsNoBg = MibiRemoveBackgroundByMaskAllChannels(countsAllSFiltCRSum,mask);

capImage=3;
for i=[3,10,11,40,42,7,12,48]
    data = countsAllSFiltCRSum(:,:,i);
    data(data>capImage)=capImage;
    figure;
    imagesc(data);
    title(massDS.Label(i));
    plotbrowser on;
end

for i=[3,10,11,40,42,7,12,48]
    data = countsNoBg(:,:,i);
    data(data>capImage)=capImage;
    figure;
    imagesc(data);
    title(massDS.Label(i));
    plotbrowser on;
end