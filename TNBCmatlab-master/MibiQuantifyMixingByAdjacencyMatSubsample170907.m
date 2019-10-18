% % MibiQuantifyMixingByAdjacencyMatSubsample170907
% Quantify using adjacent cells. For each tumor we get number of immune
% cells. We subsample each to the lowest number. We ask what percentage of
% these touch tumor cells

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
imSize=2048;
channelNum = length(massDS);
points=41;
load([pathSegment,'/dataWithTumorGroups170905.mat']);
T=50;
groupCol=55;
BootstrapNum = 100;

% % find number of immune cells per patient
% immuneNum = zeros(points,1);
% for i=1:41
%     disp(i);
%     load([pathSegment,'/Point',num2str(i),'/neighbors.mat']);
%     % reduce distance matrix only to cells classified as immune/other
%     patientInds = (dataAll(:,1) == i);
%     patientData = dataAll(patientInds,:);
%     currImmuneInds = (patientData(:,groupCol) == 2);
%     immuneNum(i) = sum(currImmuneInds);
% end
% sampleNum = min(immuneNum);

sampleNum = 131;

% for each patient subsample to sampleNum

immuneTumorInteractions = zeros(points,BootstrapNum);


for i=1:41
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/neighbors.mat']);
    % reduce distance matrix only to cells classified as immune/other
    patientInds = (dataAll(:,1) == i);
    patientData = dataAll(patientInds,:);
    currImmuneInds = (patientData(:,groupCol) == 2);
    currImuneLabels = patientData(currImmuneInds,2);
    currTumorInds = (patientData(:,groupCol) == 5) | (patientData(:,groupCol) == 6);
    currTumorLabels = patientData(currTumorInds,2);
    
    % subsample
    currImmuneNum= sum(currImmuneInds);
    parfor j=1:BootstrapNum
        randVec = randperm(currImmuneNum);
        currImuneLabelsSub = currImuneLabels(randVec([1:sampleNum]));
        immuneIndsInadjM = find(ismember([1:size(adjMat,1)],currImuneLabelsSub));
        tumorIndsInadjM = find(ismember([1:size(adjMat,1)],currTumorLabels));
        adjMatTumorImmune = adjMat(tumorIndsInadjM,immuneIndsInadjM);
        immuneTumorInteractions(i,j) = sum(sum(adjMatTumorImmune));
    end
end


z = median(immuneTumorInteractions,2);


% sort patients according to interactions percentage
[percentMixS,sInd] = sort(z);
pointsVec=[1:points];
pointsVecS = pointsVec(sInd);
figure;
bar(percentMixS);
set(gca,'XTick',[1:points],'XTickLabel',pointsVecS);

% save results
% save([pathSegment,'/mixingDataSubsample170907.mat'],'immuneTumorInteractions','pointsVecS','percentMixS');
