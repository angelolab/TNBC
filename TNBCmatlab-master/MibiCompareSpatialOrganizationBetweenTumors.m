% MibiCompareSpatialOrganizationBetweenTumors

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
imSize=2048;
channelNum = length(massDS);
points=41;
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);

%% redefine cd3-positive B cells as cd4 t-cells
%figure; histogram(dataAll(dataAll(:,57)==6,33),1000);
tcd3 = 0.5;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;

%% cell_types to compare:
% Go over all pairwise combinations of markers
%markerInds = [11,12,15:19,21:50];
markerInds = [11,12,15:19,21:28,30:50];
dataMarkers = dataAll(:,markerInds);
markerTitles = dataHeaders(markerInds);
markerNum = length(markerTitles);

zAll = zeros(markerNum,markerNum,points);

t = 1;
for i=[1:points]
    load([pathSegment,'/Point',num2str(i),'/spatialAnalysisControlTumorImmune.mat']);
    %load([pathSegment,'/Point',num2str(i),'/spatialAnalysis.mat']);
    zAll (:,:,i) = z;
end

% for all patients get the z of Ac with Tumor/Immune
cd45ind = 29;
panckind = 31;
checkInd = 25;
dataToPlot([1:points],1) = reshape(zAll(checkInd,cd45ind,:),points,1,1);
dataToPlot([1:points],2) = reshape(zAll(checkInd,panckind,:),points,1,1);

figure;
bar(dataToPlot);
title(['Compare ', markerTitles{checkInd}, ' to ', markerTitles{cd45ind}, ' and ', markerTitles{panckind}]);
legend(markerTitles{cd45ind},markerTitles{panckind});

% turn into a matrix of pairwise_comparisons x patients
spatialAll = zeros(markerNum*markerNum/2 , points);
count=0;
markerTitlesAll = {};
for i = 1 : markerNum
    for j = 1:i
        count = count+1;
        markerTitlesAll{count} = [markerTitles{i} , ' x ' , markerTitles{j}];
        for p=1:points
            spatialAll(count,p) = zAll(i,j,p);
        end
    end
end

spatialAll(isnan(spatialAll)) = 0;
spatialAll(isinf(spatialAll)) = 0;
spatialAllPlot = -spatialAll;
% remove patient 30 from the analysis
spatialAllPlotForCluster = spatialAllPlot;
spatialAllPlotForCluster(:,30) =[];
cl = clustergram(spatialAllPlotForCluster, 'ColumnLabels', [1:29,31:points], ...
        'Colormap', 'redbluecmap','DisplayRange', 0, 'DisplayRatio', 0.1);

% plot a heatmap to be able to get labels
columnLabelsN = cell2mat(cellfun(@str2num,cl.ColumnLabels,'un',0));
pairwiseIndsOrder = cell2mat(cellfun(@str2num,cl.RowLabels,'un',0));
markerTitlesAllSort = markerTitlesAll(pairwiseIndsOrder)';
figure;
imagesc(spatialAllPlot(pairwiseIndsOrder,columnLabelsN));
colormap('redbluecmap');
caxis([-20 20]);
set(gca,'xTick',[1:points],'xTickLabel',columnLabelsN,'yTick',[1:length(markerTitlesAll)], ...
    'yTickLabel',markerTitlesAllSort);
    
% plot a comparison of one marker in all patients
%
%pairst = 'Pan-Keratin x SMA';
%pairst = 'CD45 x SMA';
%pairst = 'CD45 x CD31';
%pairst = 'CD20 x CD56';
%pairst = 'CD8 x CD4';
%pairst = 'CD45 x Ki67';
pairst = {'CD45 x p53','Pan-Keratin x CD20', 'MPO x Pan-Keratin', 'PD1 x Lag3', ...
    'PD-L1 x PD1', 'IDO x PD-L1', 'CD3 x PD1', 'IDO x PD1','IDO x FoxP3'};
 
figure;
for i=1:length(pairst)
    [~ , compInd] = ismember(pairst{i},markerTitlesAll)
    vec = spatialAllPlot(compInd,columnLabelsN);
    subplot(length(pairst),1,i)
    imagesc(vec);
    colormap('redbluecmap');
    caxis([-50 50]);
    if i==length(pairst)
        set(gca,'xTick',[1:points],'xTickLabel',columnLabelsN,'yTick',1,'yTickLabel',pairst{i});
    else
        set(gca,'xTick',[],'yTick',1,'yTickLabel',pairst{i});
    end
    set(gca,'TickLength',[0 0]);
end

% plot all pairwise interactions of pd-1 across patients
marker = 'PD1';
markerAllCompIndex = find(contains(markerTitlesAll,marker));
vec = spatialAllPlot(markerAllCompIndex,columnLabelsN);
clOneMarker = clustergram(vec, 'ColumnLabels', columnLabelsN, ...
    'RowLabels',markerTitlesAll(markerAllCompIndex), ...
    'Colormap', 'redbluecmap','DisplayRange', 20, 'DisplayRatio', 0.1);

figure;
imagesc(vec);
colormap('redbluecmap');
caxis([-20 20]);
set(gca,'xTick',[1:points],'xTickLabel',columnLabelsN,'yTick',[1:size(vec,1)],'yTickLabel',markerTitlesAll(markerAllCompIndex));
set(gca,'TickLength',[0 0]);
    
    
% plot one patient
zPlot = zAll(:,:,17);
zPlot(isnan(zPlot)) = 0;
zPlot(isinf(zPlot)) = 0;
zPlot = -zPlot;
clustergram(zPlot,'RowLabels', markerTitles, 'ColumnLabels', markerTitles, ...
    'Colormap', 'redbluecmap','DisplayRange', 20, 'DisplayRatio', 0.1);


% tsne patients
spatialAllPlot2 = spatialAllPlot';
Y = tsne(spatialAllPlot2,'NumPCAComponents',29);
scatter(Y(:,1),Y(:,2));
