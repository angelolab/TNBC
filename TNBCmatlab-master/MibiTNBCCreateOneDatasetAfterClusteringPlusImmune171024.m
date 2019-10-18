% MibiTNBCCreateOneDatasetAfterClusteringPlusImmune171024.
% Add immune information on top of 

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
%path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResultsDeep';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentDeep';
immuneGroupsFile = '/Users/lkeren/Box Sync/Leeat_Share/TNBC/MapImmunePopulations171024.csv';
imSize=2048;
channelNum = length(massDS);
points=43;
categoriesImmune={'Tregs','CD4 T','CD8 T','CD3 T (gd)','NK','B','Neutrophiles','Macrophages','DC','DC/Mono','Mono/Neu','Immune'};

% add flowSOM data to immune mat
phenoS = importdata([pathCluster,'/clusterImmune171024.csv']);
flowSOMpath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/dataPerCell171019/flowsomresultsImmune171024.csv';
flowSOMS= importdata(flowSOMpath);
flowSOMVec = flowSOMS.data(:,1);
phenoS.colheaders{end+1} = 'flowSOM1';
phenoS.textdata{end+1} = 'flowSOM1';
phenoS.data(:,end+1) = flowSOMVec;

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
phenoS.data(:,1) = newPointVec;

% add immune data to dataset of all cells
load([pathSegment,'/dataWithTumorGroups171023.mat']);
groupMap = importdata(immuneGroupsFile);
% cols to add:
% 56. Immune cluster
% 57. Immune category
immuneClusterColInAll = 56;
immuneGroupColInAll = 57;

dataAll(:,immuneClusterColInAll) = 0;
dataAll(:,immuneGroupColInAll) = 0;
dataHeaders{immuneClusterColInAll} = 'immuneCluster';
dataHeaders{immuneGroupColInAll} = 'immuneGroup';

for i=1:points
    currPatientData = dataAll((dataAll(:,1) == i),:);
    currImmuneData = phenoS.data((phenoS.data(:,1) == i),:);
    immuneIndsInPdata = find(ismember(currPatientData(:,2),currImmuneData(:,2)));
    % update cluster
    currPatientData(immuneIndsInPdata,immuneClusterColInAll) = currImmuneData(:,end);
    % update category
    for k=1:length(categoriesImmune)
        clusterIDs = groupMap.data(groupMap.data(:,2) == k,1);
        for j=1:length(clusterIDs)
            clusterInds = (currPatientData(:,immuneClusterColInAll) == clusterIDs(j));
            currPatientData(clusterInds,immuneGroupColInAll) = k;
        end
    end
    % update dataAll
    dataAll((dataAll(:,1) == i),:) = currPatientData;
end

%% Fixes on top of clustering:
% 1. redefine cd3-positive B cells as cd4 t-cells
tcd3 = 0.5;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;

% 2. redefine cd68-positive Mono/Neu as Macrophages
tcd68=0.3;
dataAll( ((dataAll(:,57) == 11) & (dataAll(:,30) > tcd68)),57) = 8;

% 3. In patient 16, clusters 83 and 86 should be B cells
dataAll( ((dataAll(:,1) == 16) & ((dataAll(:,56) == 83) | (dataAll(:,56) == 86))) ,57) = 6;


save([pathSegment,'/dataWithTumorAndImmuneGroups171024.mat'],'dataAll','labelIdentityAll','dataHeaders','labelIdentityHeaders','categories','categoriesImmune');

