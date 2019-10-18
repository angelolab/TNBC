% MibiQuantifyImmunecellsPerPatientInTNBCCohort
% For all patients, get the amount of immune cells from different types and
% their distribution

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
% Plot heatmap of cells in one of the patients
markerInds = [10:12,15:30,32:50];
markerLabels = dataHeaders(markerInds);
pNum = 4;
datapNum = dataAll((dataAll(:,1) == pNum),:);
datapNum = datapNum(:,markerInds);
% CGobj = clustergram(datapNum, 'ColumnLabels', markerLabels, 'Colormap', 'redbluecmap',...
%     'Standardize',1,'RowPDist','correlation','ColumnPDist','correlation','DisplayRange',3);

% plot heatmap of adjacency matrix
load([pathSegment,'/Point',num2str(pNum),'/cellDistances.mat']);
distancesMatSmall = distancesMat([10:1110],[10:1110]);
a=colormap('Hot');
% CGobj = clustergram(distancesMatSmall, 'Colormap', a,'Symmetric','false','DisplayRange',max(max(distancesMatSmall)));

%%
% Plot clusters
clusterCol = 56;
clusterNum = max(dataAll(:,clusterCol));

% limit analysis to immune data
dataImmune = dataAll(dataAll(:,55) == 2,:);

% get median for all markers in immune clusters
markerInds = [10:12,15:30,32:50];
clusterNum = max(dataImmune(:,clusterCol));
% for each marker get the median expression in each cluster
A=dataImmune(:,clusterCol);
B=dataImmune(:,markerInds);
[xx, yy] = ndgrid(A,1:size(B,2));
markerMatMed=accumarray([xx(:) yy(:)],B(:),[],@mean);
markerLabels = dataHeaders(markerInds);

% cluster
clusterInds = [4,6:8,15,16,19,20,21,26,29,30,36];
clusterMat = markerMatMed(:,clusterInds);
clabels = cellfun(@num2str,num2cell(1:clusterNum),'uniformoutput',0);
CGobj = clustergram(clusterMat,'RowLabels', clabels , 'ColumnLabels', markerLabels(clusterInds), 'Colormap', 'redbluecmap',...
    'Standardize',1,'RowPDist','correlation','ColumnPDist','correlation','DisplayRange',3);

%% plot a heatmap with all the markers ordered by the clustering order
% get indices of markers not used for clustering
allInds = [1:length(markerLabels)];
ii = ~ismember(allInds,clusterInds);
NonClusterInds = allInds(ii);

% get numerical column labels
[~ ,clusterIndsOrder] = ismember(CGobj.ColumnLabels,markerLabels)

% get numerical row labels
%rowLabelsN = cell2mat(cellfun(@str2num,CGobj.RowLabels,'un',0));

% order clusters by final assignment
[uniqImmuneFlowSOMClusters,ia,ic] = unique(dataImmune(:,56));
clustersFinal = dataImmune(ia,57);
[clustersFinalSort sInds] = sort(clustersFinal);
rowLabelsN = uniqImmuneFlowSOMClusters(sInds);

markerMatMedCluster = markerMatMed(rowLabelsN,[clusterIndsOrder,NonClusterInds]);

% standardize the table such that in each col the mean=0 and std=1
markerMatMedClusterZ = zscore(markerMatMedCluster);
% standardize the data in each column between 0 and 1
markerMatMedClusterRescale = zeros(size(markerMatMedCluster));
for i=1:size(markerMatMedCluster,2)
    data =  markerMatMedCluster(:,i);
    dataNew = (data - min(data)) / ( max(data) - min(data) );
    markerMatMedClusterRescale(:,i) = dataNew;
end
    
% plot
figure;
imagesc(markerMatMedClusterRescale);
set(gca,'xTick',allInds,'xTickLabel',markerLabels([clusterIndsOrder,NonClusterInds]),'XTickLabelRotation',65);
set(gca,'yTick',[1:clusterNum],'yTickLabel',clabels(rowLabelsN));
%set(gca,'yTick',[]);
a = makeColorMap([1 1 1], [0 0 1], 64);
colormap(a);
set(gca,'FontSize',10);

%% Quantify immune cells per patient

countsPerPatient = zeros(length(categoriesImmune),points);
% for each patient get the total number of immune cells and their type
for i=1:points
    currInds = (dataAll(:,1)==i);
    currPatientData = dataAll(currInds,:);
    for j=1:length(categoriesImmune)
        currCategoryInds = (currPatientData(:,groupCol)==j);
        countsPerPatient(j,i) = sum(currCategoryInds);
    end
end

countsPerPatient(countsPerPatient<50) = 0;

totalPerPatient = sum(countsPerPatient,1);

% plot total
figure;
bar(totalPerPatient);

% plot normalized percentages per population
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
    0,127,255; ...
    0,0,255];
cmap01= cmap./300;

% plot total. Color by population
figure;
bar(countsPerPatient','stacked');
colormap(cmap01);

% plot pie of total
totalImmune = sum(countsPerPatient,2);
figure;
pie(totalImmune);
colormap(cmap01);

% plot pie of patient 26
figure;
pie(countsPerPatient(:,26));
colormap(cmap01);

% plot pie for patient 16
cmapPerPatient = cmap01;
patientVec = countsPerPatient(:,28);
cmapPerPatient(patientVec == 0,:) = [];
figure;
pie(patientVec);
colormap(cmapPerPatient);

% sort by total immune cells
[sVals, sInd] = sort(totalPerPatient,'descend');
countsSortedByTotal = countsPerPatient(:,sInd);

% plot bar graph of total
figure;
bar(totalPerPatient(sInd),'FaceColor',[0.5,0.5,0.5]);
axis tight;
set(gca,'xTick',[1:points],'xTickLabel',sInd,'FontSize',8,'yTick',[0,5000]);
box off;
%ylabel('Counts');

% plot total. Color by population
totalSort = totalPerPatient(sInd);
countsSortedByTotalNorm = countsSortedByTotal ./ repmat(totalSort,length(categoriesImmune),1);

figure;
% change order to make plot nicer
%changeOrderInds = [1:6,12,11,10,9,7,8];
changeOrderInds = [8,7,9,10,11,12,1:6];
countsSortedByTotalNormChangeOrder = countsSortedByTotalNorm(changeOrderInds,:);
cmap01changeOrder = cmap01(changeOrderInds,:);

subplot(1,5,[1:2]);
barh(countsSortedByTotalNormChangeOrder','stacked');
colormap(cmap01changeOrder);
set(gca,'yTick',[1:points],'yTickLabel',sInd);
ylabel('Patients sorted by total immune infiltrate','fontweight','bold');
axis tight;
title('Composition of immune cells');

subplot(1,5,3);
barh(totalSort,'FaceColor',[0.5,0.5,0.5]);
axis tight;
set(gca,'yTick',[]);
title('Total immune cells');

subplot(1,5,4);
barh(countsSortedByTotalNorm(2,:),'FaceColor',cmap01(2,:));
axis tight;
set(gca,'yTick',[]);
title('% T helper cells');

subplot(1,5,5);
barh(countsSortedByTotalNorm(8,:),'FaceColor',cmap01(8,:));
axis tight;
set(gca,'yTick',[]);
title('% Macrophages');

% get correlation for % of cd4 t-cells and macrophages with total immune
[r1 p1] = corr(totalSort',countsSortedByTotalNorm(2,:)');
[r2 p2] = corr(totalSort',countsSortedByTotalNorm(8,:)');



