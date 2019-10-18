% % MibiQuantifyImmuneTumorMixingForCohortByAdjacencyMatAndBootstrap
% Quantify using adjacent cells. For each tumor we get all the
% interactions (defined as cells touching each other). We quantify
% tumor-immune interactions out of all the interactions. We then randomize
% the labels 1000 times, and for each we get the tumor-immune interactions.
% We rank the results and the true result. The ranking of the true result
% represents the chance of seeing this amount of mixing by chance.
 

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
imSize=2048;
channelNum = length(massDS);
points=41;
load([pathSegment,'/dataWithFlowSOMGroups170823.mat']);
T=50;
BootstrapNum = 500;

immuneTumorInteractions = zeros(points,1);
immuneTumorInteractionsB = zeros(points,BootstrapNum);

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
    immuneTumorInteractions(i) = sum(sum(adjMatCellsTumorImmune));
    
end

percentMix = immuneTumorInteractions./totalInteractions;

% sort patients according to interactions percentage
[percentMixS,sInd] = sort(percentMix);
pointsVec=[1:points];
pointsVecS = pointsVec(sInd);
figure;
bar(percentMixS);
set(gca,'XTick',[1:points],'XTickLabel',pointsVecS);
