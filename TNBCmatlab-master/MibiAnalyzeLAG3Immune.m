% MibiAnalyzePD1
% for each tumor identify whether PD1 is expressed on tumor or immune
% analyze coexpression of PD1 with other markers

points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults/';
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);
groupCol = 57;

%% redefine cd3-positive B cells as cd4 t-cells
%figure; histogram(dataAll(dataAll(:,57)==6,33),1000);
tcd3 = 0.5;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;

pdl1col = find(ismember(dataHeaders,'Lag3'));

% histogram pdl1
% figure;
% histogram(dataAll(:,pdl1col),10000);

% add positive/negative pdl1 col to dataAll
t= 1;
pdl1PosVec = zeros(length(dataAll),1);
pdl1PosVec(dataAll(:,pdl1col)>t) = 1;
dataAll(:,58) = pdl1PosVec;
dataHeaders{58} = 'pd1-pos';

% For each tumor collect:
% # of immune cells
% # of pd1-positive cells
pHeaders = {'Number of Immune cells','Immune cold','Number of pdl1-pos cells','Number of pdl1-pos tumor cells','Number of pdl1-pos immune cells'};
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
    % get # of pdl1-pos cells
    pData(i,3) = sum(currData(:,58));
    % get # of pdl1-pos tumor cells
    pData(i,4) = sum( currData(:,58) & ( (currData(:,55) == 5) | (currData(:,55) == 6)));
    % get # of pdl1-pos immune cells
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

pHeaders{7} = 'Type (cold/compartmentalized/mixed)';

% remove patient 30 from the analysis
dataAll(dataAll(:,1) == 30,:) =[];
pData(30,:) = 0;

% plot # of tumor and immune pdl1 expression
% sort by tumor pdl1 expression
[sPdl1Tumor sIndsPdl1Tumor] = sort(pData(:,5),'descend');
currColor = colormap();
figure;
b=barh([pData(sIndsPdl1Tumor,5),pData(sIndsPdl1Tumor,4)],1);
b(1).FaceColor = [1 0.5 0];
b(2).FaceColor = [0 0 1];
set(gca,'yTick',[1:points],'yTickLabel',sIndsPdl1Tumor);
ylabel('Patients, sorted by PD1 positive immune cells');
xlabel('Number of PD1 positive cells');
legend({'Immune','Tumor'});
axis tight;
box off;

% plot ratio of tumor and Immune PD1
tumorPD1total = sum(pData(:,4));
immunePD1total = sum(pData(:,5));
immuneTotal = sum(pData(:,1));
tumorTotal = sum((dataAll(:,55) == 5) | (dataAll(:,55) == 6));
percentPosImmune = immunePD1total/immuneTotal;
percentPosTumor = tumorPD1total/tumorTotal;
bar([percentPosImmune,percentPosTumor],'FaceColor',[0.5 0.5 0.5]);
set(gca,'xTick',[1:2],'xTickLabel',{'Immune','Tumor'});
box off;
axis tight;

% plot ratio of Immune to tumor pd1 cells
% limit analysis to patients which have over 100 pd1 positive cells
totalPos = pData(:,4) + pData(:,5);
tumorPD1Up = pData(:,4);
tumorPD1Up(tumorPD1Up < 1) = 1;
ratioTI = tumorPD1Up ./ pData(:,5);
ratioTI(totalPos < 100) = 0;
ratioTI(isnan(ratioTI)) = 0;
ratioTI(isinf(ratioTI)) = 0;
[sRatioTI sIndsRatioTI] = sort(ratioTI,'descend');
figure;
bar(log2(sRatioTI),'FaceColor',[0.5,0.5,0.5]);
set(gca,'xTick',[1:points],'xTickLabel',sIndsRatioTI);
xlabel('Patients');
ylabel({'Ratio of PD-L1 positive cells','[log_2(tumor/immune)]'});
axis tight;
box off;

%% get data only for immune cells

cmap=[127,0,255; ...
    255,0,255; ...
    255,0,127; ...
    255,0,0; ...
    255,127,0; ...
    255,255,0; ...
    127,255,0; ...
    0,255,0; ...
    0,255,127; ...
    0,255,255; ...
    0,127,255; ];
%    0,0,255];
cmap01= cmap./300;


dataImmune = dataAll(dataAll(:,55) == 2,:);

% plot immune clusters for pdl-1 positive immune cells
pdl1PosImmune = dataImmune(dataImmune(:,pdl1col)>t,:);
pdl1NegImmune = dataImmune(dataImmune(:,pdl1col)<=t,:);

% for each patient get the composition of pd1-pos cells
pdl1PosImmuneClustersPerPatient = zeros(points,length(categoriesImmune));
for i=[1:29,31:points]
    patientPdl1PosImmune = pdl1PosImmune((pdl1PosImmune(:,1) == i),:);
    for j=1:length(categoriesImmune)
        pdl1PosImmuneClustersPerPatient(i,j) = sum((patientPdl1PosImmune(:,57) == j));
    end
end

figure;
bar(pdl1PosImmuneClustersPerPatient,'stacked');
colormap(cmap01);

% plot only t-cells
pdl1PosTcellsPerPatient = pdl1PosImmuneClustersPerPatient(:,[1:4]);
totalPdl1PosTcellsPerPatient = sum(pdl1PosTcellsPerPatient,2);
[sTotalPdl1PosTcellsPerPatient,sIndsTotalT] = sort(totalPdl1PosTcellsPerPatient,'ascend');
figure;
cmap01trunc = cmap01([1:4],:);
colormap(cmap01trunc);
b=barh(pdl1PosTcellsPerPatient(sIndsTotalT,:),'stacked');
set(gca,'yTick',[1:points],'yTickLabels',sIndsTotalT);
ylabel('Patients');
xlabel('PD1 positive T cells');
axis tight;
box off;

%% plot ratio between CD4 and CD8 cells
% limit the analysis to patients that have over 10 cells
over10Inds = find(totalPdl1PosTcellsPerPatient > 10);
cd4PosNum = sum(pdl1PosImmuneClustersPerPatient(:,[1:2]),2);
cd8PosNum = pdl1PosImmuneClustersPerPatient(:,3);
cd4PosNum(cd4PosNum == 0) = 1;
cd8PosNum(cd8PosNum == 0) = 1;
ratioCD4CD8pos = cd4PosNum ./ cd8PosNum;
ratioCD4CD8posFilt = ratioCD4CD8pos(over10Inds);
[sRatioCD4CD8pos , ratioCD4CD8posInds] = sort(ratioCD4CD8posFilt, 'descend');
type = pData(:,7);
typeOver10 = type(over10Inds)
sType = typeOver10(ratioCD4CD8posInds);

figure;
xVec = [1:length(sRatioCD4CD8pos)];
hold on;
coldInds = (sType == 2);
bar(xVec(coldInds),log2(sRatioCD4CD8pos(coldInds)),'FaceColor',[0 0 1]);
compInds = (sType == 1);
bar(xVec(compInds),log2(sRatioCD4CD8pos(compInds)),'FaceColor',[1 0 0]);
infInds = (sType == 0);
bar(xVec(infInds),log2(sRatioCD4CD8pos(infInds)),'FaceColor',[1 0.5 0]);
set(gca,'xTick',[1:length(sRatioCD4CD8pos)],'xTickLabels',over10Inds(ratioCD4CD8posInds));
xlabel('Patients');
ylabel('PD1 positive CD4 cells / PD1 positive CD8 cells');
axis tight;
box off;

figure;
b= bar(log2(sRatioCD4CD8pos));
set(gca,'xTick',[1:length(sRatioCD4CD8pos)],'xTickLabels',over10Inds(ratioCD4CD8posInds));
xlabel('Patients');
ylabel('PD1 positive CD4 cells / PD1 positive CD8 cells');
axis tight;
box off;

% check significance by ranksum test
X= sRatioCD4CD8pos(compInds);
Y= sRatioCD4CD8pos(infInds);
[p h] = ranksum(X,Y)

%% get pd1 positive b cells in specific patients
i=17;
patientIPdl1PosB = pdl1PosImmune((pdl1PosImmune(:,1) == i) & (pdl1PosImmune(:,57) == 6),:);
figure; 
subplot(3,1,1);
hpos=histogram(patientIPdl1PosB(:,56),[1:101]);
title(num2str(i));
patientIPdl1NegB = pdl1NegImmune((pdl1NegImmune(:,1) == i) & (pdl1NegImmune(:,57) == 6),:);
subplot(3,1,2);
hneg=histogram(patientIPdl1NegB(:,56),[1:101]);
ratio=hpos.BinCounts./hneg.BinCounts;
subplot(3,1,3);
bar(ratio);


% get clusters for all patients together
clusterPos = splitapply(@sum,ones(length(pdl1PosImmune),1),pdl1PosImmune(:,57));
clusterNeg = splitapply(@sum,ones(length(pdl1NegImmune),1),pdl1NegImmune(:,57));
clusterPosNorm = clusterPos/sum(clusterPos);
clusterNegNorm = clusterNeg/sum(clusterNeg);

%% plot classes of PDL1 cells. Color by population
figure;
bar([clusterPosNorm,clusterNegNorm]','stacked');
colormap(cmap01);
set(gca,'xTickLabel',{'PD1 positive','PD1 negative'});

% plot the fraction of PD1-positive cells
% partition immune cells to subtypes and also plot tumor cells
figure;
PosFraction = clusterPos./(clusterPos+clusterNeg);
PosFractionAll = [PosFraction ; percentPosTumor];
categoriesAll = [categoriesImmune([1:end-1]) , 'Other Immune' , 'Tumor'];
[sVals sInds] = sort(PosFractionAll,'descend');
bar(sVals, 'FaceColor', [0.5 0.5 0.5]);
set(gca,'xTickLabel',categoriesAll(sInds),'XTickLabelRotation',20);
ylabel('Fraction of PD1 positive cells');
axis tight;
box off;

% plot the ratio sort
figure;
PosNegRatio = clusterPosNorm./clusterNegNorm;
[sVals sInds] = sort(PosNegRatio,'descend');
bar(log2(sVals), 'FaceColor', [0.5 0.5 0.5]);
set(gca,'xTickLabel',categoriesImmune(sInds),'XTickLabelRotation',20);
ylabel('log2(PD1 positive/PD1 negative');
axis tight;
box off;

% do chi2 test to quantify pairwise stength of association
pAllchi = zeros(length(categoriesImmune),1);
for i=1:length(categoriesImmune)
    n1 = clusterPos(i);
    N1 = sum(clusterPos);
    n2 = clusterNeg(i);
    N2 = sum(clusterNeg); 
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2)
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    % Chi-square test, by hand
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    chi2stat = sum((observed-expected).^2 ./ expected);
    p = 1 - chi2cdf(chi2stat,1);
    pAllchi(i) = p;
end

figure;
boxplot(pdl1PosImmune(:,21),pdl1PosImmune(:,57));

figure;
colormap(cmap01);
pie(clusterPosNorm);
figure;
colormap(cmap01);
pie(clusterNegNorm);
legend(categoriesImmune);

%% Analyze clusters
figure;
id = 16;
pdl1id = 12;
scatter(markerMatMed(:,id),markerMatMed(:,pdl1id),'jitter','on', 'jitterAmount',0.2);
xlabel(markerLabels{id});
ylabel(markerLabels{pdl1id});

%% Scatter markers on all cells
figure;
scatCol1 = 17;
scatCol2 = 16;
dataX = dataImmune(:,scatCol1);
dataY = dataImmune(:,scatCol2);
jitterAmount = 0.2;
jitterValuesX = 2*(rand(size(dataX))-0.5)*jitterAmount;   % +/-jitterAmount max
jitterValuesY = 2*(rand(size(dataY))-0.5)*jitterAmount;   % +/-jitterAmount max
scatter(dataX+jitterValuesX, dataY+jitterValuesY);
xlabel(dataHeaders{scatCol1});
ylabel(dataHeaders{scatCol2});

%% get pdl-1 positive immune cells in patient 5
pNum =5;
pdl1PosPatient = pdl1PosImmune((pdl1PosImmune(:,1)==pNum),:);
figure;
histogram(pdl1PosPatient(:,57),12);
% get the clusters of the pdl1-pos cd8 cells
popInd=3;
pdl1PosPatientPop = pdl1PosPatient(pdl1PosPatient(:,57)==popInd,:);
figure;
histogram(pdl1PosPatientPop(:,56),100);

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
ylabel('PD-L1 positive tumor cells');
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
ylabel('PD-L1 positive immune cells');
% test significance
[pInfCold,hInfCold] = ranksum(pData((pData(:,7) == 0),4),pData((pData(:,7) == 2),4))
[pInfComp,hInfComp] = ranksum(pData((pData(:,7) == 0),4),pData((pData(:,7) == 1),4))
[pCompCold,hCompCold] = ranksum(pData((pData(:,7) == 1),4),pData((pData(:,7) == 2),4))


% plot # of tumor and immune pdl1 expression
figure;
bar([pData(:,4),pData(:,5)]);

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
