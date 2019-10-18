% MibiAnalyzePDL1
% for each tumor identify whether PDL1 is expressed on tumor or immune
% analyze coexpression of PDL1 with other markers

points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults/';
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);
groupCol = 57;

%% redefine cd3-positive B cells as cd4 t-cells
tcd3 = 0.5;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;
%%

pdl1col = find(ismember(dataHeaders,'Ki67'));

% histogram pdl1
% figure;
% histogram(dataAll(:,pdl1col),10000);

% add positive/negative pdl1 col to dataAll
%t= -0.28;
t = 1;
pdl1PosVec = zeros(length(dataAll),1);
pdl1PosVec(dataAll(:,pdl1col)>t) = 1;
dataAll(:,58) = pdl1PosVec;
dataHeaders{58} = 'Ki67-pos';

% For each tumor collect:
% # of immune cells
% # of ido-positive cells
pHeaders = {'Number of Immune cells','Immune cold','Number of ido-pos cells','Number of ido-pos tumor cells','Number of ido-pos immune cells'};
pData = zeros(points,length(pHeaders));
tImmune = 500;
for i=1:points
    pInds = (dataAll(:,1) == i);
    currData = dataAll(pInds,:);
    % get immune cells number
    immuneInds = (currData(:,55) == 2);
    pData(i,1) = sum(immuneInds);
    % mark patients with low levels of immune cells as cold
    if pData(i,1) > tImmune
        pData(i,2) = 1;
    end
    % get # of ido-pos cells
    pData(i,3) = sum(currData(:,58));
    % get # of ido-pos tumor cells
    pData(i,4) = sum( currData(:,58) & ( (currData(:,55) == 5) | (currData(:,55) == 6)));
    % get # of ido-pos immune cells
    pData(i,5) = sum( currData(:,58) & (currData(:,55) == 2));
end

% add mixing data
load([pathSegment,'/mixingDataDivByImmune170907.mat']);
%load([pathSegment,'/mixingDataSubsample170907.mat']);
[~, sInds] = sort(pointsVecS);
pData(:,6) = percentMixS(sInds);
pHeaders{6} = 'Percent mix';
% divide to compartmentalized/mixed
pData(:,7) = 0; % mixed
pData(pointsVecS([1:15]),7) = 1; %coompartmentalized
pData(pData(:,2) == 0,7) = 2; %cold
% change the classification of patient 30 to mixed. This is an error due to
% the overlap of the markers
pData(30,7) = 0; % mixed

pHeaders{7} = 'Type (cold/compartmentalized/mixed)';

%% plots
% plot # of immune cells per tumor. Color cold, infiltrating and
% compartmentalized
[sImmuneNum sIndsImmuneNum] = sort(pData(:,1),'descend');
figure;
hold on;
barh(sImmuneNum);
set(gca,'yTick',[1:41],'yTickLabel',sIndsImmuneNum);
axis tight;
xlabel('Number of immune cells');
ylabel('Patients, sorted by number of immune cells');
% cold
coldInds = find((pData(:,7) == 2));
plotColdInds = find(ismember(sIndsImmuneNum,coldInds));
barh(plotColdInds,sImmuneNum(plotColdInds),'b');
% compartmentalized
compInds = find((pData(:,7) == 1));
plotCompInds = find(ismember(sIndsImmuneNum,compInds));
barh(plotCompInds,sImmuneNum(plotCompInds),'r');
% infiltrating
infInds = find((pData(:,7) == 0));
plotInfInds = find(ismember(sIndsImmuneNum,infInds));
barh(plotInfInds,sImmuneNum(plotInfInds),'FaceColor',[1 .5 0]);

% plot box plot and data for tumor cells
figure;
hold on;
boxplot(pData(:,4),pData(:,7));
for i=1:3
    hold on;
    inds = (pData(:,7) == i-1);
    xVec = repmat(i,length(find(inds)),1);
    scatter(xVec,pData(inds,4),'filled','jitter','on', 'jitterAmount',0.2,'MarkerFaceColor',[.5 .5 .5]);
end
ylabel('IDO positive tumor cells');
% test significance
[pInfCold,hInfCold] = ranksum(pData((pData(:,7) == 0),4),pData((pData(:,7) == 2),4))
[pInfComp,hInfComp] = ranksum(pData((pData(:,7) == 0),4),pData((pData(:,7) == 1),4))
[pCompCold,hCompCold] = ranksum(pData((pData(:,7) == 1),4),pData((pData(:,7) == 2),4))

% plot box plot and data for immune cells
figure;
hold on;
boxplot(pData(:,5),pData(:,7));
for i=1:3
    hold on;
    inds = (pData(:,7) == i-1);
    xVec = repmat(i,length(find(inds)),1);
    scatter(xVec,pData(inds,5),'filled','jitter','on', 'jitterAmount',0.2,'MarkerFaceColor',[.5 .5 .5]);
end
ylabel('IDO positive immune cells');
% test significance
[pInfCold,hInfCold] = ranksum(pData((pData(:,7) == 0),4),pData((pData(:,7) == 2),4))
[pInfComp,hInfComp] = ranksum(pData((pData(:,7) == 0),4),pData((pData(:,7) == 1),4))
[pCompCold,hCompCold] = ranksum(pData((pData(:,7) == 1),4),pData((pData(:,7) == 2),4))


% plot # of tumor and immune ido expression
% sort by tumor ido expression
[sPdl1Tumor sIndsPdl1Tumor] = sort(pData(:,4),'descend');
currColor = colormap();
clr = [0 0 1;
   1 0.5 0];
figure;
colormap(clr);
barh([pData(sIndsPdl1Tumor,4),pData(sIndsPdl1Tumor,5)],1);
set(gca,'yTick',[1:points],'yTickLabel',sIndsPdl1Tumor);
ylabel('Patients, sorted by PD-L1 positive tumor cells');
xlabel('Number of IDO positive cells');
legend({'Tumor','Immune'});
axis tight;
box off;
colormap(currColor);

% plot ratio of tumor to immune pd-l1 cells
ratioTI = pData(:,4) ./ pData(:,5);
ratioTI(isnan(ratioTI)) = 0;
ratioTI(isinf(ratioTI)) = 0;
[sRatioTI sIndsRatioTI] = sort(ratioTI,'descend');
figure;
bar(log2(sRatioTI),'FaceColor',[0.5,0.5,0.5]);
set(gca,'xTick',[1:points],'xTickLabel',sIndsRatioTI);
xlabel('Patients');
ylabel({'Ratio of IDO positive cells','[log_2(tumor/immune)]'});
axis tight;
box off;

% sort by tumor pdl1 expression
% [sPdl1Tumor sIndsPdl1Tumor] = sort(pData(:,4));
% figure;
% hold on;
% bar(log10(sPdl1Tumor));
% set(gca,'xTick',[1:41],'xTickLabel',sIndsPdl1Tumor);
% % plot cold tumors
% coldInds = find(pData(:,2) == 0);
% [~,coldsIndsS] = ismember(sIndsPdl1Tumor,coldInds);
% bar(find(coldsIndsS),log10(sPdl1Tumor(find(coldsIndsS))),'r');

% same plot. Mark by different colors the X coldest tumors/
% compartmentalized tumors / mixed tumors
numColor = 7;
[sPdl1Tumor sIndsPdl1Tumor] = sort(pData(:,4));
figure;
hold on;
bar(log10(sPdl1Tumor),0.8);
set(gca,'xTick',[1:41],'xTickLabel',sIndsPdl1Tumor);
% plot cold tumors
[~, sColdInds] = sort(pData(:,1));
sColdIndsToTake = sColdInds([1:numColor]);
[~,coldIndsInS] = ismember(sIndsPdl1Tumor,sColdIndsToTake);
bar(find(coldIndsInS),log10(sPdl1Tumor(find(coldIndsInS))),0.8,'EdgeColor','r','FaceColor','r');
% % plot compartmentalized
% % first remove negative inds
% nonNegativeInds = [1:points];
% nonNegativeInds(sColdIndsToTake) =[];
% [~, sCompInds] = sort(pData(nonNegativeInds,6));
% sCompIndsToTake = sCompInds([1:numColor]);
% pointsToTake = nonNegativeInds(sCompIndsToTake);
% [~,compIndsInS] = ismember(sIndsPdl1Tumor,pointsToTake);
% % get bar width
% plotXInds = find(compIndsInS);
% divVal = min(plotXInds([2:end])-plotXInds([1:end-1]));
% bw = 0.8/divVal;
% bar(find(compIndsInS),log10(sPdl1Tumor(find(compIndsInS))),bw,'EdgeColor','g','FaceColor','g');
% plot mixed
[~, sMixInds] = sort(pData(nonNegativeInds,6),'descend');
sMixIndsToTake = sMixInds([1:numColor]);
pointsToTake = nonNegativeInds(sMixIndsToTake);
[~,mixIndsInS] = ismember(sIndsPdl1Tumor,pointsToTake);
% get bar width
plotXInds = find(mixIndsInS);
divVal = min(plotXInds([2:end])-plotXInds([1:end-1]));
bw = 0.8/divVal;
bar(find(mixIndsInS),log10(sPdl1Tumor(find(mixIndsInS))),bw,'EdgeColor','c','FaceColor','c');
xlabel('Patients (Ranked by number of PDL1-positive tumor cells)');
ylabel('log10[Number of PDL1-positive tumor cells');

% do a ranksum test to gauge significance
pCold = ranksum(sPdl1Tumor(find(coldIndsInS)),sPdl1Tumor(find(~coldIndsInS)));
pMix = ranksum(sPdl1Tumor(find(mixIndsInS)),sPdl1Tumor(find(~mixIndsInS)));

%% plot # of immune cells vs. mixing
figure;
scatter(pData(:,1), pData(:,6))
xlabel('# of Immune cells');
ylabel('mixing coefficient');

%% sort patients by # of immune cells, and plot pdl1 tumor cells
figure;
[sImmuneNum sIndsImmuneNum] = sort(pData(:,1));
figure;
bar(pData(sIndsImmuneNum,4));
set(gca,'xTick',[1:41],'xTickLabel',sIndsImmuneNum);

%% sort patients by # of immune cells, and plot pdl1 immune cells
figure;
[sImmuneNum sIndsImmuneNum] = sort(pData(:,1));
figure;
bar(pData(sIndsImmuneNum,5));
set(gca,'xTick',[1:41],'xTickLabel',sIndsImmuneNum);
