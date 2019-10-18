points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/';
pathCluster = [path,'ClusteringResults/'];
resultsDir = [path,'dataPerCell170801/Cluster1/'];

phenoS = importdata([pathCluster,'/170803_exportPhenoDataStd.csv']);
phenoCol = 58;
tsne1Col = 53;
tsne2Col = 54;
clusterInds = [12,15:19,22,24,26:27,30:33,35:43,45,48:49];
clusterNum = max(phenoS.data(:,phenoCol));
% for each marker get the median expression in each cluster
A=phenoS.data(:,phenoCol);
B=phenoS.data(:,clusterInds);
[xx, yy] = ndgrid(A,1:size(B,2));
clusterMat=accumarray([xx(:) yy(:)],B(:),[],@median);

% cluster
clabels = cellfun(@num2str,num2cell(1:clusterNum),'uniformoutput',0);
CGobj = clustergram(clusterMat,'RowLabels', clabels , 'ColumnLabels', phenoS.textdata(clusterInds), 'Colormap', 'redbluecmap',...
    'Standardize',1,'RowPDist','correlation','ColumnPDist','correlation','DisplayRange',3);

% classify clusters to categories
cellCategories = {'Immune','SMA','Tumor','Endothel'};
%categoryInds = {[25,16,11,30,10,9,12,13,20,29,6,23,18,26],[14,8,24],[2,5,3,22,31,28,19,21,27,4,1,7,17],[15]};
categoryInds = {[16,11,9,12,13,36,33,7,35,20,34,24,31,21],[29,3,28,15],[25,19,32,5,2,22,27,8,30,14,18,26,4,1,23,10,6],[17]};

categoryVec= zeros(length(phenoS.data),1);
for i=1:length(cellCategories)
    categoryClusters = categoryInds{i};
    for j=1:length(categoryClusters)
        categoryVec(phenoS.data(:,phenoCol) == categoryClusters(j)) = i;
    end
end

dataMat = [phenoS.data(:,[2:end]) , categoryVec];
Labels2Write = [ phenoS.colheaders([2:end]) , 'cellCategory' ];

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

% for each cell add a category label and save fcs files and data structures
% for all cells, for immune and for tumor
mkdir(resultsDir);
TEXT.PnS = Labels2Write;
TEXT.PnN = Labels2Write;

for i=1:points
    dataPoint = dataMat((newPointVec == i),:);
    dataImmune = dataPoint((dataPoint(:,end) == 1),:);
    dataTumor = dataPoint((dataPoint(:,end) == 3),:);
    save([path,'/Point',num2str(i),'/cellDataCategory.mat'],'dataPoint','Labels2Write');
    writeFCS([resultsDir,'/dataStdAll_p',num2str(i),'.fcs'],dataPoint,TEXT); 
    writeFCS([resultsDir,'/dataStdImmune_p',num2str(i),'.fcs'],dataImmune,TEXT);
    writeFCS([resultsDir,'/dataStdTumor_p',num2str(i),'.fcs'],dataTumor,TEXT);
end
