% MibiTNBC plot heatmap of flowSOM results
points=43;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResultsDeep';
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
phenoS = importdata([pathCluster,'/clusterImmune171024.csv']);
phenoCol = 55;
resultsDir = [path,'/dataPerCell171019/ClusterImmune171024/'];


% add flowSOM data
flowSOMpath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/dataPerCell171019/flowsomresultsImmune171024.csv';
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
set(gca,'FontSize',10);
a = makeColorMap([1 1 1], [0 0 1], 64);
colormap(a);

% plot a bar graph with the percentage of cells in each cluster
clusterPercentages = histcounts(phenoS.data(:,phenoCol),[1:101])./length(phenoS.data);
figure;
bar(clusterPercentages);

% % look at pdl1 cd11c correlation
% figure;
% id = 6;
% pdl1id = 12;
% scatter(markerMatMed(:,id),markerMatMed(:,pdl1id),'jitter','on', 'jitterAmount',0.2);
% xlabel(markerLabels{id});
% ylabel(markerLabels{pdl1id});
% 
% markerLabelsSort = markerLabels([clusterIndsOrder,NonClusterInds]);
% id=9;
% figure; scatter(markerMatMedClusterRescale(:,id),markerMatMedClusterRescale(:,21));
% [r p]=corr(markerMatMedClusterRescale(:,id),markerMatMedClusterRescale(:,21));
% markerLabelsSort=markerLabels([clusterIndsOrder,NonClusterInds]);
% xlabel(markerLabelsSort(id));
% ylabel(markerLabelsSort(21));

% % all cells
% id=18;
% figure;
% scatplot(phenoS.data(:,id), phenoS.data(:,23));
% %scatter(phenoS.data(:,id),phenoS.data(:,23));
% [r p]=corr(phenoS.data(:,id),phenoS.data(:,23));
% xlabel(phenoS.colheaders(id));
% ylabel('pdl1');


%% analyze specific population
popNum=83;

% get the patient distribution for this cluster

% convert values in 1st colums to point number
% cyt saves the order sorted as a string, so the numbers don't represent
% our points
pointStr = cellfun(@num2str,num2cell([1:29,31:points+1]),'uniformoutput',0);
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
% immuneCategories = [23,34,33,47,36,25,14,4,13,44,46,35,24,22,31,21,11,1,2,3,12,32,41,42,45,16,6,26,5,15,17,7,74,56,37,18,8,9,19,27,38,28,10,20, ...
%     48,55,29,49,39,30,40,50];
% allCategories = [1:clusterNum];
% otherCategories = setdiff(allCategories,immuneCategories);
% categoryInds = {immuneCategories,otherCategories};
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

% % fixes
% 1. fix classification for cluster 41 in patient 30 from immune to tumor
% %%%%%%%%%%%%%%%
% fixPatient = 30;
% fixCluster = 41;
% fixVal = 2;
% fixInds = ((newPointVec == fixPatient) & (phenoS.data(:,phenoCol)==fixCluster));
% dataMat(fixInds,[end]) = fixVal;
% 
% 2. fix cluster 58, which was 10% of the data classified as negative
% %%%%%%%%%%%%%%%%
% fixCluster = 58;
% badInds = (phenoS.data(:,phenoCol)==fixCluster);
% badData = dataMat(badInds,:);
% plot cd45 vs. keratin vs. beta for this cluster
% m=[badData(:,42),badData(:,44)]';
% [d, h] = densityScatter(m);
% figure;
% scatter(badData(:,42),badData(:,44),10,newPointVec(badInds));
% xlabel('cd45');
% ylabel('panCK');
% cells which have cd45 higher than T45 and panCK lower than TpanCK should be classified accordingly
% T45 = -0.5;
% TpanCK = -0.4;
% cd45Col = 43;
% panCKcol = 45;
% fixVal = 1;
% fixInds = ((phenoS.data(:,phenoCol)==fixCluster) & (phenoS.data(:,cd45Col)>T45) & (phenoS.data(:,panCKcol)<TpanCK));
% dataMat(fixInds,[end]) = fixVal; 


% %for each cell add a category label and save fcs files and data structures
% %for all cells, for immune and for tumor
% mkdir(resultsDir);
% TEXT.PnS = Labels2Write;
% TEXT.PnN = Labels2Write;
% 
% for i=1:points
%     dataPoint = dataMat((newPointVec == i),:);
%     dataImmune = dataPoint((dataPoint(:,end) == 1),:);
%     dataTumor = dataPoint((dataPoint(:,end) == 2),:);
%     save([path,'/Point',num2str(i),'/cellDataCategory170831.mat'],'dataPoint','Labels2Write');
%     writeFCS([resultsDir,'/dataStdAll_p',num2str(i),'.fcs'],dataPoint,TEXT); 
%     writeFCS([resultsDir,'/dataStdImmune_p',num2str(i),'.fcs'],dataImmune,TEXT);
%     writeFCS([resultsDir,'/dataStdTumor_p',num2str(i),'.fcs'],dataTumor,TEXT);
% end
% 
% % generate a struct for all the data to save
% phenoSOld = phenoS;
% phenoS.data = [newPointVec , dataMat];
% phenoS.textdata = ['GateSource' , Labels2Write];
% phenoS.colheaders = ['GateSource' , Labels2Write];
% 
% %save([pathCluster,'/dataSeparateImmuneTumor170831_3.mat'],'phenoS');