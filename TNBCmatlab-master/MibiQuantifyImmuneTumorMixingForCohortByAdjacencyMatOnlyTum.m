% % MibiQuantifyImmuneTumorMixingForCohortByAdjacencyMat
% Quantify using adjacent cells. For each tumor we get all the
% interactions (defined as cells touching each other). We quantify
% tumor-immune interactions out of all the interactions and that is the
% mixing percentage
% problem with this approach - it is sensitive to the number of
% immune/tumor cells in the tissue
% 

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
imSize=2048;
channelNum = length(massDS);
points=41;
load([pathSegment,'/dataWithFlowSOMGroups170823.mat']);
T=50;

immuneTumorInteractions = zeros(points,1);
tumorNum = zeros(points,1);

for i=1:points
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/cellDistances.mat']);
    load([pathSegment,'/Point',num2str(i),'/neighbors.mat']);
    % reduce distance matrix only to cells classified as immune/other
    currInds = (labelIdentityAll(:,1) == i);
    currLabelIdentity = labelIdentityAll(currInds,:);
    cellInds = currLabelIdentity((currLabelIdentity(:,3)==1),2);
    adjMatCells = adjMat(cellInds,cellInds);
    % Get interactions between Tumor and immune cells
    currInds2 = (dataAll(:,1) == i);
    currCellData = dataAll(currInds2,:);
    % get tumor and immune inds
    currTumorInds = currCellData(:,54) == 0;
    currImmuneInds = currCellData(:,54) > 0;
    adjMatCellsTumorImmune = adjMatCells(currTumorInds>0,currImmuneInds>0);
    % get number of Tumor cells
    tumorNum(i) = length(find(currTumorInds));
    % get number of tumor cells touching immune cells
    sumTumor = sum(adjMatCellsTumorImmune,2);
    sumTumor(sumTumor>1)=1; % this step is to move from interactions to number of positive cells
    immuneTumorInteractions(i) = sum(sumTumor);
end

percentMix = immuneTumorInteractions./tumorNum;

% sort patients according to interactions percentage
[percentMixS,sInd] = sort(percentMix);
pointsVec=[1:points];
pointsVecS = pointsVec(sInd);
figure;
bar(percentMixS);
set(gca,'XTick',[1:points],'XTickLabel',pointsVecS);
