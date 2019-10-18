% MibiTNBCCreateOneDatasetAfterClustering 170901.
% Add information after separating tumor and immune cells, gating tumor and
% clustering it

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
%path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResultsDeep';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentDeep';
restGroupsFile = '/Users/lkeren/Box Sync/Leeat_Share/TNBC/MapTumorRestPopulations171023.csv';
imSize=2048;
channelNum = length(massDS);
points=43;
categories={'Unidentified','Negative','Immune','Endothel','Mesenchimal-like','Tumor','Keratin-pos Tumor'};

load([pathCluster,'/dataGateTumor171023.mat']);
phenoSTumorRest = importdata([pathCluster,'/clusterTumorRest171023.csv']);
groupMap = importdata(restGroupsFile);
% tumor column inds
phenoS_TumorImmuneCol = 54;
% tumor rest column inds
phenoSTumorRest_ClusterCol = 56;
cellIdCol = 2;

dataAll = zeros(1,channelNum+6);
labelIdentityAll = zeros(1,1);

% convert values in 1st colums to point number
% cyt saves the order sorted as a string, so the numbers don't represent
% our points
pointStr = cellfun(@num2str,num2cell([1:29,31:points+1]),'uniformoutput',0);
pointStrSort=sort(pointStr);
pointsByStrCell = cellfun(@str2num,pointStrSort,'uniformoutput',0);
pointsByStr = cell2mat(pointsByStrCell);

newPointVec = zeros(length(phenoSTumorRest.data),1);
for i=1:points
    pointFromCyt = i;
    realPoint = pointsByStr(i);
    newPointVec(phenoSTumorRest.data(:,1) == pointFromCyt) = realPoint;
end
phenoSTumorRest.data(:,1) = newPointVec;

% cols in the final dataset:
% 1. Patient
% 2. Cell ID
% 3. Size
% 4-53. data
% 54. TumorYN
% 55. Tumor cluster
% 56. Combined category: {'Unidentified','Negative','Immune','Endothel','Mesenchimal-like','Tumor','Keratin-pos Tumor'}

% unite the clustering data with the tumor data
ind=1;
indIdentity=1;
for i=[1:29,31:points+1]
    load([pathSegment,'/Point',num2str(i),'/cellData.mat']);
    load([pathSegment,'/Standard171019/Point',num2str(i),'/cellData.mat']);
    dataCurr=[];
    currNumCells = length(labelVec);
    dataCurr([1:currNumCells],1) = i;
    dataCurr([1:currNumCells],2) = labelVec;
    dataCurr([1:currNumCells],3) = cellSizesVec;
    dataCurr([1:currNumCells],4:channelNum+4-1) = dataScaleSizeCellsTransStd;
    dataCurr([1:currNumCells],channelNum+4) = 0;
    dataCurr([1:currNumCells],channelNum+5) = 0;
    dataCurr([1:currNumCells],channelNum+6) = 0;
    % add tumor-immune information
    patientTumorInds = (phenoS.data(:,1) == i);
    patientTumorData = phenoS.data(patientTumorInds,:);
    cellLabels = patientTumorData(:,cellIdCol);
    insertInds = find(ismember(dataCurr(:,2),cellLabels));
    dataCurr(insertInds,channelNum+4) = 1; % this column marks that it was initially classified as not immune
    % add tumor Rest cluster information
    patientTumorRestInds = (phenoSTumorRest.data(:,1) == i);
    patientTumorRestData = phenoSTumorRest.data(patientTumorRestInds,:);
    cellLabelsRest = patientTumorRestData(:,cellIdCol);
    insertIndsRest = find(ismember(dataCurr(:,2),cellLabelsRest));
    dataCurr(insertIndsRest,channelNum+5) = patientTumorRestData(:,phenoSTumorRest_ClusterCol); % this column has the tumor rest cluster
    % add final information
    % 1 - mark negative cells
    clusterIDs = groupMap.data(groupMap.data(:,2) == 1,1);
    for j=1:length(clusterIDs)
        clusterInds = (dataCurr(:,channelNum+5) == clusterIDs(j));
        dataCurr(clusterInds,channelNum+6) = 1;
    end
    % 2 - mark immune cells
    immuneInds = (dataCurr(:,channelNum+4) == 0);
    dataCurr(immuneInds,channelNum+6) = 2;
    % 3 - mark cd31 cells
    clusterIDs = groupMap.data(groupMap.data(:,2) == 2,1);
    for j=1:length(clusterIDs)
        clusterInds = (dataCurr(:,channelNum+5) == clusterIDs(j));
        dataCurr(clusterInds,channelNum+6) = 3;
    end
    % 4 - mark mesenchimal-like cells
    clusterIDs = groupMap.data(groupMap.data(:,2) == 3,1);
    for j=1:length(clusterIDs)
        clusterInds = (dataCurr(:,channelNum+5) == clusterIDs(j));
        dataCurr(clusterInds,channelNum+6) = 4;
    end
    % 5 - mark beta-catenin tumor
    clusterIDs = groupMap.data(groupMap.data(:,2) == 4,1);
    for j=1:length(clusterIDs)
        clusterInds = (dataCurr(:,channelNum+5) == clusterIDs(j));
        dataCurr(clusterInds,channelNum+6) = 5;
    end
    % 6 - mark tumor keratin cells
    keratinInds = ((dataCurr(:,channelNum+4) == 1) & (dataCurr(:,channelNum+5) == 0));
    dataCurr(keratinInds,channelNum+6) = 6;

    % add patient to big table
    dataAll([ind:ind+currNumCells-1],:) = dataCurr;
    ind = ind+currNumCells;
    labelIdentityAll([indIdentity:indIdentity+length(labelIdentityNew2)-1],1) = i;
    labelIdentityAll([indIdentity:indIdentity+length(labelIdentityNew2)-1],2) = [1:length(labelIdentityNew2)];
    labelIdentityAll([indIdentity:indIdentity+length(labelIdentityNew2)-1],3) = labelIdentityNew2;
    indIdentity = indIdentity+length(labelIdentityNew2);
end

dataHeaders = ['SampleID'; channelLabelsForFCS; 'tumorYN'; 'tumorCluster' ; 'Group'];
labelIdentityHeaders = {'SampleID';'cellLabelInImage';'identity'};

save([pathSegment,'/dataWithTumorGroups171023.mat'],'dataAll','labelIdentityAll','dataHeaders','labelIdentityHeaders','categories');

