% MibiPlotImmunePopulationsAndHighlightsForFigure2
% plot classification of cells to immune groups for all patients

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults/';
imSize=2048;
channelNum = length(massDS);
points=41;
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);
T=50;
tumorYNcol = 53;
groupCol=57;

%% redefine cd3-positive B cells as cd4 t-cells
tcd3 = 0;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;

images=zeros(imSize,imSize,3,points,'uint8');

% colorcode:
% [ 0 bg - white, 0 tumor  - gray, immune- rainbow ]
cmap = [200 200 200; 127 0 255; 255 0 255; 255 0 127; 255 0 0; 255 127 0; 255 255 0; 127 255 0; 0 255 0; 0 255 127; 0 255 255; 0 127 255; 0 0 255];
cmap01 = cmap/255;

rangeX = [1070:1350];
rangeY = [250:530]; 

for i=[12];
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/segmentationParams.mat']);
    stats = regionprops(newLmod,'PixelIdxList');
    currInds = (dataAll(:,1) == i);
    currCellData = dataAll(currInds,:);
    currIndsLabel = (labelIdentityAll(:,1) == i);
    currLabelIdentity = labelIdentityAll(currIndsLabel,:);
    labelNum = length(currLabelIdentity);
    imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
    counts=0;
    for j=1:labelNum-1
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % check if tumor
            if (currCellData(cellInd,tumorYNcol) == 1)
                cellVal =1;
            else
                cellVal = currCellData(cellInd,groupCol) + 1;
            end
            imageL(stats(j).PixelIdxList)=cellVal;
        end
    end
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    images(:,:,:,i) = label2rgb(imageLReshape,cmap01,'w');
    % plot single
    f=figure;
    colormap('parula');
    imagesc(label2rgb(imageLReshape,cmap01,'w'));
    axis square;
    % plot subset
    f2 = figure;
    imagesc(label2rgb(imageLReshape(rangeX,rangeY),cmap01,'w'));
    axis square;
    set(gca,'xTick',[],'yTick',[]);
%    plotbrowser on;
%     saveas(f,[path,'/Point',num2str(i),'/TIFs3/immuneGroups.png']);
end

% plot individual channels
% for i=[12];
%     disp(i);
%     channelsToPlot = [12,14,29,35,16,27];
%     load([path,'/Point',num2str(i),'/dataDeNoiseCohortNoAgg.mat']);
%     for j = channelsToPlot
%         currMat = countsNoNoiseNoAgg(rangeX,rangeY,j);
%         capVal = 3;
%         currMat(currMat>capVal) = capVal;
%         currMatFilt = imgaussfilt(currMat,1);
%         figure;
%         colormap('gray');
%         imagesc(currMatFilt);
%         axis square;
%         set(gca,'xTick',[],'yTick',[]);
%     end
% end

% plot color overlays
channelsToOverlay = [12,14];
for i=[12];
    matToPlot=zeros(length(rangeX),length(rangeY),3);
    for j = 1:length(channelsToOverlay)
        currMat = countsNoNoiseNoAgg(rangeX,rangeY,channelsToOverlay(j));
        capVal = 3;
        currMat(currMat>capVal) = capVal;
        currMatFilt = imgaussfilt(currMat,1);
        matToPlot(:,:,j) = currMatFilt;
    end
    figure;
    imshow(matToPlot);
end


% get order to plot
% load([pathSegment,'/mixingDataDivByTotal170907.mat']);
% 
% plot all together
% f=figure;
% [ha, pos] = tight_subplot(6,7,[.01 .01],[.01 .01],[.01 .01]); 
% for j = 1:points
%     axes(ha(j));
%     data = images(:,:,:,j);
%     imshow(data);
%     title(['Sample ',num2str(j)]);
% end
% saveas(f,[path,'/ByChannel3/groupsByImmuneMixing170911.fig']);