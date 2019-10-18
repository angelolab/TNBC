% MibiPrepareMultiLayerTiffForHaloForCohort

points=41;
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
mkdir([path,'/HALO/']);
processNucleiFlag =1;

for pointNumber = 2:5
    load([path,'/Point',num2str(pointNumber),'/dataDeNoiseCohortNoAgg.mat']);
    % add a nuclear channel summimg the nuclear markers
    channelLabels = massDS.Label;
    channelLabels{end+1} = 'nuclei';
    channels = {'dsDNA','H3K27me3','H3K9ac'};
    [tf chInd] = ismember(channels,massDS.Label);
    nucData = countsNoNoiseNoAgg(:,:,chInd);
    nucDataSum = sum(nucData,3);
    % save images of nuclei
    if processNucleiFlag
        nucDataGaus = imgaussfilt(nucDataSum,1);
        nucDataGray = mat2gray(nucDataGaus);
        nucDataAdj = imadjust(nucDataGray);
        imwrite(nucDataAdj,[path,'/HALO/NucleiGaus_p',num2str(pointNumber),'.tif']);
    else
        nucDataGray = mat2gray(nucDataSum);
        imwrite(nucDataGray,[path,'/HALO/Nuclei_p',num2str(pointNumber),'.tif']);
    
    % add nuclear channel to multi-layer tif and print
%     allData = cat(3,countsNoNoiseNoAgg,nucDataSum);
%     MibiPrepareMultiLayerTiffForHalo (allData,channelLabels,[path,'/HALO/Multi_p',num2str(pointNumber),'.tif']);
%     MibiPrepareMultiLayerTiffForHalo (countsNoNoiseNoAgg,massDS.Label,[path,'/HALO/p_',num2str(pointNumber),'_multi.tif']);
    end
end