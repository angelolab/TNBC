% MibiAnalyzeSpatialLocationOfPDL1Immune
% for each patient identify whether PDL1 immune cells are interacting with
% the tumor more than random

points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults/';
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);
imSize=2048;
channelNum = length(massDS);
groupCol = 55;
BootstrapNum = 100;
pdl1col = find(ismember(dataHeaders,'PD-L1'));

% histogram pdl1
% figure;
% histogram(dataAll(:,pdl1col),10000);

% add positive/negative pdl1 col to dataAll
t= -0.2;
pdl1PosVec = zeros(length(dataAll),1);
pdl1PosVec(dataAll(:,pdl1col)>t) = 1;
dataAll(:,58) = pdl1PosVec;
dataHeaders{58} = 'pdl1-pos';


%PDL1posImmuneTumorInteractions = zeros(points,1);
%PDL1posImmuneTumorInteractionsB = zeros(points,BootstrapNum);

for i=1:points
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/neighbors.mat']);
    % reduce distance matrix only to cells classified as immune/other
    patientInds = (dataAll(:,1) == i);
    patientData = dataAll(patientInds,:);
    currImmunePdl1PosInds = ( (patientData(:,groupCol) == 2) & (patientData(:,58) == 1) );
    currImunePdl1PosLabels = patientData(currImmunePdl1PosInds,2);
    currTumorInds = (patientData(:,groupCol) == 5) | (patientData(:,groupCol) == 6);
    currTumorLabels = patientData(currTumorInds,2);
    
    immunePdl1PosIndsInadjM = find(ismember([1:size(adjMat,1)],currImunePdl1PosLabels));
    tumorIndsInadjM = find(ismember([1:size(adjMat,1)],currTumorLabels));
    adjMatTumorImmune = adjMat(tumorIndsInadjM,immunePdl1PosIndsInadjM);
    interactingCells=sum(adjMatTumorImmune,1);
    PDL1posImmuneTumorInteractions(i) = sum(interactingCells>0);
    
    pdl1posImmuneNum = length(immunePdl1PosIndsInadjM);
    currImmuneInds = (patientData(:,groupCol) == 2);
    currImmuneLabels = patientData(currImmuneInds,2);
    immuneNum = length(currImmuneLabels);
    
    %% plot and color pdl1-pos immune cells 
%     load([pathSegment,'/Point',num2str(i),'/segmentationParams.mat']);
%     stats = regionprops(newLmod,'PixelIdxList');
%     currIndsLabel = (labelIdentityAll(:,1) == i);
%     currLabelIdentity = labelIdentityAll(currIndsLabel,:);
%     labelNum = length(currLabelIdentity);
%     imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
%     counts=0;
%     for j=1:labelNum-1
%         cellInd = find(ismember(patientData(:,2),j));
%         if ~isempty(cellInd)
%             if ismember(j,currImunePdl1PosLabels)
%                 cellVal = 1;
%             elseif ismember(j,currImmuneLabels)
%                 cellVal = 2;
%             elseif ismember(j,currTumorLabels)
%                 cellVal = 3;
%             else
%                 cellVal = 4;
%             end
%             imageL(stats(j).PixelIdxList)=cellVal;
%         end
%     end
%     imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
%     f=figure;
%     colormap('parula');
%     imagesc(label2rgb(imageLReshape));
%     plotbrowser on;
    
    %% randomize location of pdl1-pos immune cells
    for j = 1:BootstrapNum
        %randPDL1
        randPdl1PosLabels = datasample(currImmuneLabels,pdl1posImmuneNum,'Replace',false);
        randImmunePdl1PosIndsInadjM = find(ismember([1:size(adjMat,1)],randPdl1PosLabels));
        randAdjMatCellsTumorImmune = adjMat(tumorIndsInadjM,randImmunePdl1PosIndsInadjM);
        randInteractingCells=sum(randAdjMatCellsTumorImmune,1);
        PDL1posImmuneTumorInteractionsB(i,j) = sum(randInteractingCells>0);
        
%         % plot
%         imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
%         counts=0;
%         for k=1:labelNum-1
%             cellInd = find(ismember(patientData(:,2),k));
%             if ~isempty(cellInd)
%                 if ismember(k,randPdl1PosLabels)
%                     cellVal = 1;
%                 elseif ismember(k,currImmuneLabels)
%                     cellVal = 2;
%                 elseif ismember(k,currTumorLabels)
%                     cellVal = 3;
%                 else
%                     cellVal = 4;
%                 end
%                 imageL(stats(k).PixelIdxList)=cellVal;
%             end
%         end
%         imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
%         f=figure;
%         colormap('parula');
%         imagesc(label2rgb(imageLReshape));
%         plotbrowser on;
        
    end
end

var = PDL1posImmuneTumorInteractionsB;
realVar = PDL1posImmuneTumorInteractions;
for i=1:points
    [muhat(i),sigmahat(i)] = normfit(var(i,:));
    z(i) = (muhat(i)-realVar(i))/sigmahat(i);
end

% sort patients according to interactions percentage
[percentMixS,sInd] = sort(z);
pointsVec=[1:points];
pointsVecS = pointsVec(sInd);
figure;
bar(percentMixS);
set(gca,'XTick',[1:points],'XTickLabel',pointsVecS);