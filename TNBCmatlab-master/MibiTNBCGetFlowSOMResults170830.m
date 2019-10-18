% plot heatmap of flowSOM results
points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
phenoS = importdata([path,'/170830_exportImmunePhenoDataStd.csv']);
phenoCol = 69;

% add flowSOM data
flowSOMpath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/dataPerCell170809/ClusterAll170829/dataStdImmune/flowsomresults.csv';
flowSOMS= importdata(flowSOMpath);
flowSOMVec = flowSOMS.data(:,1);
phenoS.colheaders{end+1} = 'flowSOM1';
phenoS.textdata{end+1} = 'flowSOM1';
phenoS.data(:,end+1) = flowSOMVec;

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
clusterInds = [4,6:8,15,16,19,20,21,26,29,30,36];
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
clusterPercentages = histcounts(phenoS.data(:,phenoCol),[1:101])./length(phenoS.data);
figure;
bar(clusterPercentages);

%% analyze specific population
popNum=57;

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
