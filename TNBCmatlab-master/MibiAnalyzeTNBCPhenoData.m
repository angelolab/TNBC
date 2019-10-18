points=41;
% massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
% path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
% phenoS = importdata([path,'/170804_exportImmunePhenoDataStd.csv']);
% phenoCol = 69;
% tsne1Col = 67;
% tsne2Col = 68;

% massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
% path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
% phenoS = importdata([path,'/170809_exportPhenoDataStd.csv']);
% phenoCol = 63;
% tsne1Col = 57;
% tsne2Col = 58;

%% Immune cells clustered with cd11c and cd11b 170809
% massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
% path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
% phenoS = importdata([path,'/170809_exportImmunePhenoDataStd.csv']);
% phenoCol = 72;
% tsne1Col = 66;
% tsne2Col = 67;
% clusterNum = max(phenoS.data(:,phenoCol));

%% Immune cells clustered without cd11c and cd11b 170816
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
phenoS = importdata([path,'/170809_exportImmunePhenoDataStd.csv']);
phenoCol = 88;
tsne1Col = 74;
tsne2Col = 75;
clusterNum = max(phenoS.data(:,phenoCol));

% figure;
% histogram(phenoS.data(:,58));

figure;
colormap('jet');
specificCluster=ones(length(phenoS.data),1);
specificCluster((phenoS.data(:,phenoCol)) == 1) = 2;
scatter(phenoS.data(:,tsne1Col),phenoS.data(:,tsne2Col),1,specificCluster);

% for each marker get the median expression in each cluster
clusterInds = [15:19,26,27,30,32:33,38,42,48];
A=phenoS.data(:,phenoCol);
B=phenoS.data(:,clusterInds);
[xx, yy] = ndgrid(A,1:size(B,2));
clusterMat=accumarray([xx(:) yy(:)],B(:),[],@median);

% cluster
clabels = cellfun(@num2str,num2cell(1:clusterNum),'uniformoutput',0);
CGobj = clustergram(clusterMat,'RowLabels', clabels , 'ColumnLabels', phenoS.textdata(clusterInds), 'Colormap', 'redbluecmap',...
    'Standardize',1,'RowPDist','cosine','ColumnPDist','cosine','DisplayRange',3);
