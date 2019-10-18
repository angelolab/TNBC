% MibiPlotTumorImmunecellsForCohort

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
imSize=2048;
channelNum = length(massDS);
points=41;
load([pathSegment,'/dataWithFlowSOMGroups170823.mat']);
T=50;

images=zeros(imSize,imSize,3,points);

for i=1:points
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/segmentationParams.mat']);
    stats = regionprops(newLmod,'PixelIdxList');
    currInds2 = (dataAll(:,1) == i);
    currCellData = dataAll(currInds2,:);
    % get tumor and immune inds
    currTumorInds = currCellData((currCellData(:,54) == 0),2);
    currImmuneInds = currCellData((currCellData(:,54) > 0),2);
    % generate the R channel - tumor cells
    rChannel = zeros(imSize,imSize);
    % generate the G channel - imuune cells
    gChannel = zeros(imSize,imSize);
    for j=1:length(stats)
        if ismember(j, currTumorInds)
            rChannel(stats(j).PixelIdxList)=1;
        elseif ismember(j, currImmuneInds)
            gChannel(stats(j).PixelIdxList)=1;
        end 
    end
    images(:,:,1,i) = rChannel;
    images(:,:,2,i) = gChannel;
    f=figure;
    imagesc(cat(3,rChannel,gChannel,zeros(imSize)));
    title(['Point ',num2str(i)]);
    plotbrowser on;
    %saveas(f,[path,'/point',num2str(i),'/TIFs3/tumorImmune.png']);
    imwrite(cat(3,rChannel,gChannel,zeros(imSize)),[path,'/point',num2str(i),'/TIFs3/tumorImmune.tif']);
end

% plot all together
f=figure;
[ha, pos] = tight_subplot(6,7,[.01 .01],[.01 .01],[.01 .01]); 
for j = 1:points
    axes(ha(j));
    data = images(:,:,:,j);
    imshow(data);
    title(['Sample ',num2str(j)]);
end
saveas(f,[path,'/ByChannel3/tumorImmune170829.fig']);