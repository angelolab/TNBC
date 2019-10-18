% MibiTNBCMarkTumorClusters170828
% for all cells add a column marking their 'primary' population
points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/';
pathCluster = [path,'ClusteringResults/'];
resultsDir = [path,'dataPerCell170809/ClusterTumor/'];

phenoS = importdata([pathCluster,'/170901_TumorRest.csv']);
% phenoCol = 80;
% tsne1Col = 74;
% tsne2Col = 75;
phenoCol = 74;
tsne1Col = 72;
tsne2Col = 73;
clusterNum = max(phenoS.data(:,phenoCol));

% get median for all markers in clusters
markerInds = [10:12,15:30,32:50];
clusterNum = max(phenoS.data(:,phenoCol));
% for each marker get the median expression in each cluster
A=phenoS.data(:,phenoCol);
B=phenoS.data(:,markerInds);
[xx, yy] = ndgrid(A,1:size(B,2));
markerMatMed=accumarray([xx(:) yy(:)],B(:),[],@mean);
markerLabels = phenoS.textdata(markerInds);

% cluster
clusterInds = [2,3,11,28];
clusterMat = markerMatMed(:,clusterInds);
clabels = cellfun(@num2str,num2cell(1:clusterNum),'uniformoutput',0);
CGobj = clustergram(clusterMat,'RowLabels', clabels , 'ColumnLabels', markerLabels(clusterInds), 'Colormap', 'redbluecmap',...
    'Standardize',1,'RowPDist','correlation','ColumnPDist','correlation','DisplayRange',3);

%% plot a heatmap with all the markers ordered by the clustering order
% get indices of markers not used for clustering
allInds = [1:length(markerLabels)];
ii = ~ismember(allInds,clusterInds);
NonClusterInds = allInds(ii);

% get numerical column labels
[~ ,clusterIndsOrder] = ismember(CGobj.ColumnLabels,markerLabels)

% get numerical row labels
rowLabelsN = cell2mat(cellfun(@str2num,CGobj.RowLabels,'un',0));

markerMatMedCluster = markerMatMed(rowLabelsN,[clusterIndsOrder,NonClusterInds]);

% standardize the table such that in each col the mean=0 and std=1
markerMatMedClusterZ = zscore(markerMatMedCluster);
% standardize the data in each column between 0 and 1
markerMatMedClusterRescale = zeros(size(markerMatMedCluster));
for i=1:size(markerMatMedCluster,2)
    data =  markerMatMedCluster(:,i);
    dataNew = (data - min(data)) / ( max(data) - min(data) );
    markerMatMedClusterRescale(:,i) = dataNew;
end
    
% plot
figure;
imagesc(markerMatMedClusterRescale);
set(gca,'xTick',allInds,'xTickLabel',markerLabels([clusterIndsOrder,NonClusterInds]),'XTickLabelRotation',45);
set(gca,'yTick',[1:clusterNum],'yTickLabel',clabels(rowLabelsN));
a = makeColorMap([1 1 1], [0 0 1], 64);
colormap(a);

% plot a bar graph with the percentage of cells in each cluster
clusterPercentages = histcounts(phenoS.data(:,phenoCol),[1:clusterNum+1])./length(phenoS.data);
figure;
bar(clusterPercentages);

%% analyze specific population
popNum=5;
% view population on tsne

figure;
colormap('jet');
specificCluster=ones(length(phenoS.data),1);
specificCluster((phenoS.data(:,phenoCol)) == popNum) = 2;
scatter(phenoS.data(:,tsne1Col),phenoS.data(:,tsne2Col),1,specificCluster);

% get the patient distribution for this cluster

% convert values in 1st colums to point number
% cyt saves the order sorted as a string, so the numbers don't represent
% our points
pointStr = cellfun(@num2str,num2cell(1:points),'uniformoutput',0);
pointStrSort=sort(pointStr);
pointsByStrCell = cellfun(@str2num,pointStrSort,'uniformoutput',0);
pointsByStr = cell2mat(pointsByStrCell);

newPointVec = zeros(length(phenoS.data),1);
for i=1:points
    pointFromCyt = i;
    realPoint = pointsByStr(i);
    newPointVec(phenoS.data(:,1) == pointFromCyt) = realPoint;
end

popLog = (phenoS.data(:,phenoCol) == popNum);
popPatients= newPointVec(popLog);
patientCountForPop = histcounts(popPatients,[1:points+1]);
figure;
bar([1:points],patientCountForPop);
grid on;


% %% classify clusters to categories
% cellCategories = {'Immune','Other'};
% %categoryInds = {[25,16,11,30,10,9,12,13,20,29,6,23,18,26],[14,8,24],[2,5,3,22,31,28,19,21,27,4,1,7,17],[15]};
% categoryInds = {[23,33,17,1,31,41,6,3,15,34,16,7,5,9,27,36],[46,43,28,37,38,26,4,8,14,35,11,24,18,29,13,2,22,10,45,39,30,48,21,12,20,19,32,47,40,42,25,44]};
% 
% categoryVec= zeros(length(phenoS.data),1);
% for i=1:length(cellCategories)
%     categoryClusters = categoryInds{i};
%     for j=1:length(categoryClusters)
%         categoryVec(phenoS.data(:,phenoCol) == categoryClusters(j)) = i;
%     end
% end
% 
% dataMat = [phenoS.data(:,[2:end]) , categoryVec];
% Labels2Write = [ phenoS.colheaders([2:end]) , 'cellCategory' ];
% 
% % convert values in 1st colums to point number
% % cyt saves the order sorted as a string, so the numbers don't represent
% % our points
% pointStr = cellfun(@num2str,num2cell(1:points),'uniformoutput',0);
% pointStrSort=sort(pointStr);
% pointsByStrCell = cellfun(@str2num,pointStrSort,'uniformoutput',0);
% pointsByStr = cell2mat(pointsByStrCell);
% 
% newPointVec = zeros(length(phenoS.data),1);
% for i=1:points
%     pointFromCyt = i;
%     realPoint = pointsByStr(i);
%     newPointVec(phenoS.data(:,1) == pointFromCyt) = realPoint;
% end
% 
% % for each cell add a category label and save fcs files and data structures
% % for all cells, for immune and for tumor
% % mkdir(resultsDir);
% % TEXT.PnS = Labels2Write;
% % TEXT.PnN = Labels2Write;
% % 
% % for i=1:points
% %     dataPoint = dataMat((newPointVec == i),:);
% %     dataImmune = dataPoint((dataPoint(:,end) == 1),:);
% %     dataTumor = dataPoint((dataPoint(:,end) == 2),:);
% %     save([path,'/Point',num2str(i),'/cellDataCategory.mat'],'dataPoint','Labels2Write');
% %     writeFCS([resultsDir,'/dataStdAll_p',num2str(i),'.fcs'],dataPoint,TEXT); 
% %     writeFCS([resultsDir,'/dataStdImmune_p',num2str(i),'.fcs'],dataImmune,TEXT);
% %     writeFCS([resultsDir,'/dataStdTumor_p',num2str(i),'.fcs'],dataTumor,TEXT);
% % end
