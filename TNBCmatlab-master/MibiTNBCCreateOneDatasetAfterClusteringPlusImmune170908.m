% MibiTNBCCreateOneDatasetAfterClusteringPlusImmune 170908.
% Add immune information on top of 

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
%path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
immuneGroupsFile = '/Users/lkeren/Box Sync/Leeat_Share/TNBC/MapImmunePopulations170908.csv';
imSize=2048;
channelNum = length(massDS);
points=41;
categoriesImmune={'Tregs','CD4 T','CD8 T','CD3 T (gd)','NK','B','Neutrophiles','Macrophages','DC','DS/Mono','Mono/Neu','Immune'};

% add flowSOM data to immune mat
phenoS = importdata([pathCluster,'/170907_ClusterImmune.csv']);
flowSOMpath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/dataPerCell170809/ClusterAll170831/flowsomresultsImmune170907.csv';
flowSOMS= importdata(flowSOMpath);
flowSOMVec = flowSOMS.data(:,1);
phenoS.colheaders{end+1} = 'flowSOM1';
phenoS.textdata{end+1} = 'flowSOM1';
phenoS.data(:,end+1) = flowSOMVec;

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
phenoS.data(:,1) = newPointVec;

% add immune data to dataset of all cells
load([pathSegment,'/dataWithTumorGroups170905.mat']);
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

save([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat'],'dataAll','labelIdentityAll','dataHeaders','labelIdentityHeaders','categories','categoriesImmune');

