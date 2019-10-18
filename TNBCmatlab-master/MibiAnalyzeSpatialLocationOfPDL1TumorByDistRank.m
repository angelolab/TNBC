% MibiAnalyzeSpatialLocationOfPDL1ImmuneBy Dist
% for each patient identify whether PDL1 tumor cells are interacting with
% the immune cells more than random. We rank all tumor cells by distance to
% imune and then look at enrichment of pdl1 positive cells

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

closeRatio = zeros(points,1);
closePdl1 = zeros(points,1);
closeAll = zeros(points,1);
farRatio = zeros(points,1);
farPdl1 = zeros(points,1);
farAll = zeros(points,1);
chiScore = zeros(points,1);

for i=1:points
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/cellDistances.mat']);
    % reduce distance matrix only to cells classified as immune/other
    patientInds = (dataAll(:,1) == i);
    patientData = dataAll(patientInds,:);
    currImmuneInds = (patientData(:,groupCol) == 2);
    currImmuneLabels = patientData(currImmuneInds,2);
    currTumorInds = (patientData(:,groupCol) == 5) | (patientData(:,groupCol) == 6);
    currTumorLabels = patientData(currTumorInds,2);
    currTumorData = patientData(currTumorInds,:);
    tumorNum = length(currTumorLabels);
    % get for each tumor cell the distance to the 5 closest immune cells
    adjMatTumorImmune = distancesMat(currImmuneLabels,currTumorLabels);
    k=5;
    sAdjMatTumorImmune = sort(adjMatTumorImmune);
    sAdjMatTumorImmuneK = sAdjMatTumorImmune([1:k],:);
    distToImmune = mean(sAdjMatTumorImmuneK);
    closeInds = distToImmune < 100;
    closeLabels = currTumorData(closeInds,2);
    
    % plot and color immune cells with mean distance <100
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
%             if ismember(j,closeLabels)
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
    
    % count number of pd-l1 positive cells in close and far cells
    closePos = closeInds' & (currTumorData(:,58) == 1);
    closePdl1(i) = sum(closePos);
    closeAll(i) = sum(closeInds);
    closeRatio(i) = sum(closePos) / sum(closeInds);
    farPos = ~closeInds' & (currTumorData(:,58) == 1);
    farPdl1(i) = sum(farPos);
    farAll(i) = sum(~closeInds);
    farRatio(i) = sum(farPos) / sum(~closeInds);
    
    % perfor chi2 test to assess significance
   n1 = closePdl1(i); N1 = closeAll(i);
   n2 = farPdl1(i); N2 = farAll(i);
   % Pooled estimate of proportion
   p0 = (n1+n2) / (N1+N2);
   % Expected counts under H0 (null hypothesis)
   n10 = N1 * p0;
   n20 = N2 * p0;
   % Chi-square test, by hand
   observed = [n1 N1-n1 n2 N2-n2];
   expected = [n10 N1-n10 n20 N2-n20];
   chi2stat = sum((observed-expected).^2 ./ expected);
   p = 1 - chi2cdf(chi2stat,1);
    
   chiScore(i) = p;
    
end

% limit the analysis to patients which had PD-L1 expression predominantly
% on tumor cells
pImmuneList = [33,27,8,23,31,21,26,18,25,7,14,29,38,36,20,3,2,17,11,40,12,13,16];
closeRatioFilt = closeRatio(pImmuneList);
farRatioFilt = farRatio(pImmuneList);
chiFilt = chiScore(pImmuneList);

% limit the analysis to patients with >500 PD-L1 positive cells
tpdl1 = 500;
pdl1All = closePdl1+farPdl1;
pdl1Filt = pdl1All(pImmuneList);
filt2Inds = (pdl1Filt > tpdl1);

pImmuneListFilt2 = pImmuneList(filt2Inds);
closeRatioFilt2 = closeRatioFilt(filt2Inds);
farRatioFilt2 = farRatioFilt(filt2Inds);
chiFilt2 = chiFilt(filt2Inds);

% sort by the differences in ratio
closFarFilt2ratio = closeRatioFilt2./farRatioFilt2;
[sclosFarFilt2ratio sInds] = sort(closFarFilt2ratio,'descend');
chiFilt2Sort = chiFilt2(sInds);

figure;
b=bar([closeRatioFilt2(sInds) , farRatioFilt2(sInds)]);
set(gca,'xTick',[1:length(pImmuneListFilt2)],'xTickLabel',pImmuneListFilt2(sInds));
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0.5 0.5];
xlabel('Patients');
ylabel('Ratio of PD-L1 positive cells');
legend({'Close to immune','Distant from immune'});
box off;


% var = log2(closeRatio./farRatio);
% [sVar sInds] = sort(var);
% figure;
% bar(sVar);
% set(gca,'xTick',[1:points],'xTickLabel',sInds);