% MibiUniteFlowSOMClustersTNBCCohort
% Add group annotations to FlowSOM clusters

points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
phenoS = importdata([path,'/170809_exportImmunePhenoDataStd.csv']);
phenoCol = 92;

% add flowSOM data
flowSOMpath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/dataPerCell170809/ClusterImmune/flowsomresults.csv';
flowSOMS= importdata(flowSOMpath);
flowSOMVec = flowSOMS.data(:,1);
phenoS.colheaders{end+1} = 'flowSOM1';
phenoS.textdata{end+1} = 'flowSOM1';
phenoS.data(:,end+1) = flowSOMVec;

clusterNum = max(phenoS.data(:,phenoCol));

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

phenoS.colheaders{end+1} = 'patientNum';
phenoS.textdata{end+1} = 'patientNum';
phenoS.data(:,end+1) = newPointVec;

% group clusters
groupsS = importdata('/Users/lkeren/Box Sync/Leeat_Share/TNBC/MapImmunePopulations.csv');
originalCluster = groupsS.data(:,1);
group = groupsS.data(:,2);

% generate a vector with the group per patient

newMap = containers.Map(num2cell(originalCluster),num2cell(group));
groupVec = zeros(length(newPointVec),1);
for i=1:length(newPointVec)
    groupVec(i) = newMap(phenoS.data(i,phenoCol));
end

phenoS.colheaders{end+1} = 'group1';
phenoS.textdata{end+1} = 'group1';
phenoS.data(:,end+1) = groupVec;

save([path,'/ClusterWithGroups170821.mat'],'phenoS');
