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
load([pathSegment,'/dataWithTumorGroups170905.mat']);
T=50;
groupCol=55;
BootstrapNum = 500;

immuneTumorInteractions = zeros(points,1);
totalInteractions = zeros(points,1);
immuneInteractions = zeros(points,1);
tumorInteractions = zeros(points,1);

% immuneTumorInteractionsB = zeros(points,BootstrapNum);
% immuneInteractionsB = zeros(points,BootstrapNum);
% tumorInteractionsB = zeros(points,BootstrapNum);

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
    
    immuneIndsInadjM = find(ismember([1:size(adjMat,1)],currImuneLabels));
    tumorIndsInadjM = find(ismember([1:size(adjMat,1)],currTumorLabels));
    adjMatTumorImmune = adjMat(tumorIndsInadjM,immuneIndsInadjM);
    adjMatTumor = adjMat(tumorIndsInadjM,tumorIndsInadjM);
    adjMatImmune = adjMat(immuneIndsInadjM,immuneIndsInadjM);
    adjMatTogether = adjMat([tumorIndsInadjM,immuneIndsInadjM],[tumorIndsInadjM,immuneIndsInadjM]);
    immuneTumorInteractions(i) = sum(sum(adjMatTumorImmune));
    tumorInteractions(i) = sum(sum(adjMatTumor))/2;
    immuneInteractions(i) = sum(sum(adjMatImmune))/2;
    totalInteractions(i) = sum(sum(adjMatTogether))/2;
    
%     % randomize immune and tumor
%     tumorNum = length(tumorIndsInadjM);
%     immuneNum = length(immuneIndsInadjM);
%     totalNum = tumorNum+immuneNum;
%     parfor j = 1:BootstrapNum
%         randVec = randperm(totalNum);
%         randTumorInds = randVec([1:tumorNum]);
%         randImmuneInds = randVec([tumorNum+1:end]);
%         randAdjMatCellsTumorImmune = adjMatTogether(randTumorInds,randImmuneInds);
%         randAdjMatCellsTumor = adjMatTogether(randTumorInds,randTumorInds);
%         randAdjMatCellsImmune = adjMatTogether(randImmuneInds,randImmuneInds);
%         immuneTumorInteractionsB(i,j) = sum(sum(randAdjMatCellsTumorImmune));
%         tumorInteractionsB(i,j) = sum(sum(randAdjMatCellsTumor))/2;
%         immuneInteractionsB(i,j) = sum(sum(randAdjMatCellsImmune))/2;
%     end
end

% var = immuneTumorInteractionsB./(immuneInteractionsB+immuneTumorInteractionsB);
% realVar = immuneTumorInteractions./(immuneInteractions+immuneTumorInteractions);
% for i=1:points
%     [muhat(i),sigmahat(i)] = normfit(var(i,:));
%     z(i) = (muhat(i)-realVar(i))/sigmahat(i);
% end

z = immuneTumorInteractions./(immuneInteractions+immuneTumorInteractions);
z = immuneTumorInteractions./(tumorInteractions+immuneTumorInteractions);
%z = immuneTumorInteractions./totalInteractions;

% sort patients according to interactions percentage
[percentMixS,sInd] = sort(z);
pointsVec=[1:points];
pointsVecS = pointsVec(sInd);
figure;
bar(percentMixS);
set(gca,'XTick',[1:points],'XTickLabel',pointsVecS);

% % save randomizations
% save([pathSegment,'/mixingData170907.mat'],'immuneTumorInteractionsB','immuneTumorInteractions','immuneInteractionsB','immuneInteractions','tumorInteractionsB','tumorInteractions','totalInteractions');
% 
% save results
save([pathSegment,'/mixingDataDivByTumor170907.mat'],'immuneTumorInteractions','z','pointsVecS','percentMixS');
