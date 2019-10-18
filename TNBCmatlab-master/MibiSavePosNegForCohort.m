% MibiSavePosNegForCohort
% for each cell mark whether it is pos or neg for each marker.
% for each patient get the number of positive cells for each marker

points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults/';
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);
groupCol = 57;

categoriesImmune{10} = 'DC/Mono';

%% redefine cd3-positive B cells as cd4 t-cells
%figure; histogram(dataAll(dataAll(:,57)==6,33),1000);
tcd3 = 0.5;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;

pd1col = find(ismember(dataHeaders,'PD1'));
pdl1col = find(ismember(dataHeaders,'PD-L1'));
lag3col = find(ismember(dataHeaders,'Lag3'));
ido1col = find(ismember(dataHeaders,'IDO'));

% add positive/negative column to all markers
t= 1;
posNegCellTable = zeros(length(dataAll),length(massDS));
dataAllOnlyMarkers = dataAll(:,[4:52]);
posNegCellTable(dataAllOnlyMarkers>t) = 1;

% For each tumor collect:
% # of positive cells for each marker
% # of immune cells
% # of pd-l1 positive tumor cells
% # of pd-l1 positive immune cells

pData = zeros(points,length(massDS)+7);
pHeaders = [massDS.Label ; 'Number of Immune cells' ; 'Cold'; 'Number of pd-l1-pos tumor cells' ; 'Number of pd-l1-pos immune cells'; ...
    'Ratio of pdl1 tumor and immune' ; 'Number of pd1-pos cd4 cells' ; 'Number of pd1-pos cd8 cells' ; ...
    'Ratio of pd1-pos CD4 and CD8 cells' ; 'Number of lag3-pos and pd1-pos cells' ; 'Number of lag3-pos Tregs' ; ...
    'Number of lag3-pos cd4 cells' ; 'Number of lag3-pos cd8 cells' ; 'Number of Tregs' ; 'Number of lag3-neg Tregs' ; ...
    'Number of CD4 T cells' ; 'Number of CD8 T cells' ; 'Number of pd1-pos lag3-pos cd8' ; 'Number of pd1-pos lag3-pos cd4' ; ...
    'Number of pd1-pos T cells' ; 'Number of lag3-pos T cells' ; 'Number of pd1-lag3-pos T cells' ; 'Number of T cells' ; ...
    'Number of tumor cells' ; 'Number of ido-pos cells' ; 'Number of ido-pos immune cells' ; 'Number of ido-pos tumor cells' ; ...
    'Number of pd1-pos cells' ; 'Number of pd1-pos immune cells' ; 'Number of pd1-pos tumor cells' ; ...
    'Number of lag3-pos cells' ; 'Number of lag3-pos immune cells' ; 'Number of lag3-pos tumor cells' ; ...
    'Number of ki67-pos tumor cells' ; 'Number of p53-pos tumor cells' ; 'Number of egfr-pos tumor cells' ; ...
    'Number of ck6-pos tumor cells' ; 'Number of ck17-pos tumor cells' ; 'Number of hla1-pos tumor cells' ; ...
    ' Number of hladr-pos cells']; 
tImmune = 500;
for i=1:points
    pInds = (dataAll(:,1) == i);
    currData = dataAll(pInds,:);
    currPosNeg = posNegCellTable(pInds,:);
    % get positive cells for each marker
    for j=1: length(massDS)
        markerVec = zeros(length(currData),1);
        markerVec(currData(:,j+3) > t) = 1;
        pData(i,j) = sum(markerVec);
    end
    % get immune cells number
    immuneInds = (currData(:,55) == 2);
    pData(i,length(massDS)+1) = sum(immuneInds);
    % mark patients with low levels of immune cells as cold
    if pData(i,length(massDS)+1) > tImmune
        pData(i,length(massDS)+2) = 1;
    end
    % get # of pdl1-pos tumor cells
    currPdlPosInds = (currData(:,pdl1col) > t);
    pData(i,length(massDS)+3) = sum( currPdlPosInds & ( (currData(:,55) == 5) | (currData(:,55) == 6)));
    % get # of pdl1-pos immune cells
    pData(i,length(massDS)+4) = sum( currPdlPosInds & (currData(:,55) == 2));
    % get ratio of tumor and immune
    tumorNum = pData(i,length(massDS)+3);
    tumorNum(tumorNum ==0) = 1;
    immuneNum = pData(i,length(massDS)+4);
    immuneNum(immuneNum ==0) = 1;
    pData(i,length(massDS)+5) = tumorNum/immuneNum;
    % get # of pd1-pos cd4 cells
    currPd1PosInds = (currData(:,pd1col) > t);
    pData(i,length(massDS)+6) = sum( currPd1PosInds & ( (currData(:,57) == 1) | (currData(:,57) == 2)));
    % get # of pd1-pos cd8 cells
    currPd1PosInds = (currData(:,pd1col) > t);
    pData(i,length(massDS)+7) = sum( currPd1PosInds & (currData(:,57) == 3) );
    % get ratio of PD1 positive CD8 and CD4 cells
    % make 0 values into 1 to get data
    cd4Num = pData(i,length(massDS)+6);
    cd4Num(cd4Num == 0) =1;
    cd8Num = pData(i,length(massDS)+7);
    cd8Num(cd8Num == 0) =1;
    pData(i,length(massDS)+8) = cd8Num/cd4Num;
    % get # of lag3-pos + pd1-pos cells
    currLag3PosInds = (currData(:,lag3col) > t);
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+9) = sum( currPd1PosInds & currLag3PosInds);
    end
    % get # of lag3-pos treg cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+10) = sum( currLag3PosInds & (currData(:,57) == 1));
    end
    % get # of lag3-pos cd4 cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+11) = sum( currLag3PosInds & (currData(:,57) == 2));
    end
    % get # of lag3-pos cd8 cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+12) = sum( currLag3PosInds & (currData(:,57) == 3) );
    end
    % get number of Tregs
    pData(i,length(massDS)+13) = sum(currData(:,57) == 1);
    % get number of lag3-neg tregs
    pData(i,length(massDS)+14) = pData(i,length(massDS)+13)-pData(i,length(massDS)+10);
    % get number of CD4 t cells
    pData(i,length(massDS)+15) = sum(currData(:,57) == 2);
    % get number of CD8 t cells
    pData(i,length(massDS)+16) = sum(currData(:,57) == 3);
    % get number of lag3-pos and pd1-pos cd8 t cells
    currLag3PosInds = (currData(:,lag3col) > t);
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+17) = sum( currPd1PosInds & currLag3PosInds & (currData(:,57) == 3));
    end
    % get number of lag3-pos and pd1-pos cd4 t cells
    currLag3PosInds = (currData(:,lag3col) > t);
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+18) = sum( currPd1PosInds & currLag3PosInds & (currData(:,57) == 2));
    end
    % get number of pd1-pos t cells
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currPd1PosInds)
        pData(i,length(massDS)+19) = sum( currPd1PosInds & ...
            ( (currData(:,57) == 1) | (currData(:,57) == 2) | (currData(:,57) == 3) | (currData(:,57) == 4) ) );
    end
    % get number of lag3-pos t cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+20) = sum( currLag3PosInds & ...
            ( (currData(:,57) == 1) | (currData(:,57) == 2) | (currData(:,57) == 3) | (currData(:,57) == 4) ) );
    end
    % get number of lag3-pos and pd1-pos t cells
    currLag3PosInds = (currData(:,lag3col) > t);
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+21) = sum( currPd1PosInds & currLag3PosInds & ...
            ( (currData(:,57) == 1) | (currData(:,57) == 2) | (currData(:,57) == 3) | (currData(:,57) == 4) ) );
    end
    % get number of t cells
    pData(i,length(massDS)+22) = sum( (currData(:,57) == 1) | (currData(:,57) == 2) | (currData(:,57) == 3) | (currData(:,57) == 4) );
    % get number of tumor cells
    tumorInds = (currData(:,55) == 5) | (currData(:,55) == 6);
    pData(i,length(massDS)+23) = sum(tumorInds);
    % get number of ido-pos cells
    currIdoPosInds = (currData(:,ido1col) > t);
    if ~isempty(currIdoPosInds)
        pData(i,length(massDS)+24) = sum( currIdoPosInds );
    end
    % get number of ido-pos immune cells
    currIdoPosInds = (currData(:,ido1col) > t);
    if ~isempty(currIdoPosInds)
        pData(i,length(massDS)+25) = sum( currIdoPosInds & (currData(:,55) == 2));
    end
    % get number of ido-pos tumor cells
    currIdoPosInds = (currData(:,ido1col) > t);
    if ~isempty(currIdoPosInds)
        pData(i,length(massDS)+26) = sum( currIdoPosInds & ( (currData(:,55) == 5) | (currData(:,55) == 6) ) );
    end
    % get number of pd1-pos cells
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currPd1PosInds)
        pData(i,length(massDS)+27) = sum( currPd1PosInds );
    end
    % get number of pd1-pos immune cells
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currPd1PosInds)
        pData(i,length(massDS)+28) = sum( currPd1PosInds & (currData(:,55) == 2));
    end
    % get number of pd1-pos tumor cells
    currPd1PosInds = (currData(:,pd1col) > t);
    if ~isempty(currPd1PosInds)
        pData(i,length(massDS)+29) = sum( currPd1PosInds & ( (currData(:,55) == 5) | (currData(:,55) == 6) ) );
    end
        % get number of lag3-pos cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+30) = sum( currLag3PosInds );
    end
    % get number of lag3-pos immune cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+31) = sum( currLag3PosInds & (currData(:,55) == 2));
    end
    % get number of lag3-pos tumor cells
    currLag3PosInds = (currData(:,lag3col) > t);
    if ~isempty(currLag3PosInds)
        pData(i,length(massDS)+32) = sum( currLag3PosInds & ( (currData(:,55) == 5) | (currData(:,55) == 6) ) );
    end
end

% add mixing data
load([pathSegment,'/mixingDataDivByImmune170907.mat']);
%load([pathSegment,'/mixingDataSubsample170907.mat']);
[~, sInds] = sort(pointsVecS);
pData(:,length(massDS)+33) = percentMixS(sInds);
pHeaders{length(massDS)+33} = 'Percent mix';
% divide to compartmentalized/mixed
pData(:,length(massDS)+34) = 0; % mixed
pData(pointsVecS([1:20]),length(massDS)+34) = 1; %coompartmentalized
pData(pData(:,length(massDS)+2) == 0,length(massDS)+34) = 2; %cold

pHeaders{length(massDS)+34} = 'Type (cold/compartmentalized/mixed)';

removeInds = [1:7,10,41,43,48,49];
pData(:,removeInds) = [];
pHeaders(removeInds) = [];

% Heatmap
figure;
imagesc(pData);
set(gca,'xTick',[1:length(pHeaders)],'xTickLabel',pHeaders,'XTickLabelRotation',40);
%clustergram(pData,'ColumnLabels',pHeaders,'Symmetric',false,'Standardize',1);
%indsVec = [10,12,40,41,42,43,44,45];
%indsVec = [1:38,40:length(pHeaders)];
indsVec = [1:length(pHeaders)];
indsVec=[5,40:51];
clustergram(pData([1:29,31:points],indsVec),'ColumnLabels',pHeaders(indsVec),'Symmetric',false,'Standardize',1 ...
    ,'DisplayRatio',0.1,'ColumnLabelsRotate',35,'RowLabels',[1:29,31:points],'Colormap','redbluecmap', 'RowPDist','euclidean','ColumnPDist','euclidean');

%% plot bar graphs of immune checkpoint molecules in tumor/immune
% Order: PD1 , Lag3 , PD-L1 , IDO
% Order: Immune , tumor
checkPointNames ={'PD-1','LAG3','PD-L1','IDO'};
checkPointCols = [pd1col, lag3col, pdl1col, ido1col];
immuneTumorCheckpoint = zeros(4,2);
immuneTumorCheckpoint(1,1) = sum(pData(:,65)) ./ sum(pData(:,38));
immuneTumorCheckpoint(1,2) = sum(pData(:,66)) ./ sum(pData(:,60));
immuneTumorCheckpoint(2,1) = sum(pData(:,68)) ./ sum(pData(:,38));
immuneTumorCheckpoint(2,2) = sum(pData(:,69)) ./ sum(pData(:,60));
immuneTumorCheckpoint(3,1) = sum(pData(:,41)) ./ sum(pData(:,38));
immuneTumorCheckpoint(3,2) = sum(pData(:,40)) ./ sum(pData(:,60));
immuneTumorCheckpoint(4,1) = sum(pData(:,62)) ./ sum(pData(:,38));
immuneTumorCheckpoint(4,2) = sum(pData(:,63)) ./ sum(pData(:,60));

figure;
colormap([1 0.5 0 ; 0 0 1]);
bar(immuneTumorCheckpoint);
set(gca,'xTick',[1:length(immuneTumorCheckpoint)],'xTickLabels',checkPointNames);
ylabel('Fraction of positive cells');
box off;
axis tight;

%%
% For each checkpoint molecule get expression in each one of the immune
% clusters + tumor
dataTumor = dataAll((dataAll(:,55) == 5) | (dataAll(:,55) == 6),:);
dataImmune = dataAll(dataAll(:,55) == 2,:);
cpPosImmuneClustersPerPatient = zeros(points,length(categoriesImmune),length(checkPointNames));
% go over all chackpoints
for k=1:length(checkPointNames)
    currCol = checkPointCols(k);
    cpPosImmune = dataImmune(dataImmune(:,currCol)>t,:);
    cpNegImmune = dataImmune(dataImmune(:,currCol)<=t,:);
    % for each patient get the composition of checkpoint-pos cells
    for i=[1:29,31,33:points]
        patientCpPosImmune = cpPosImmune((cpPosImmune(:,1) == i),:);
        for j=1:length(categoriesImmune)
            cpPosImmuneClustersPerPatient(i,j,k) = sum((patientCpPosImmune(:,57) == j));
        end
    end
    % get clusters for all patients together
    clusterPos(:,k) = splitapply(@sum,ones(length(cpPosImmune),1),cpPosImmune(:,57));
    clusterNeg(:,k) = splitapply(@sum,ones(length(cpNegImmune),1),cpNegImmune(:,57));
    clusterPosNorm(:,k) = clusterPos(:,k)/sum(clusterPos(:,k));
    clusterNegNorm(:,k) = clusterNeg(:,k)/sum(clusterNeg(:,k));
end

% get baseline composition for all imune clusters
dataImmuneNo_30_32 = dataImmune(((dataImmune(:,1) ~= 30) & (dataImmune(:,1) ~= 32)),:);
clusterCompAll = splitapply(@sum,ones(length(dataImmuneNo_30_32),1),dataImmuneNo_30_32(:,57));
clusterCompAllNorm = clusterCompAll / sum(clusterCompAll);
% plot total composition
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
figure;
bar([clusterCompAllNorm,clusterPosNorm]','stacked');
colormap(cmap01);
set(gca,'xTick',[1:length(checkPointNames)+1],'xTickLabel',['Baseline',checkPointNames],'XTickLabelRotation',30);
box off;

% for each checkpoint in each subtype plot ratio between percentage and
% baseline percentage
clusterPosRatio = clusterPosNorm ./ repmat(clusterCompAllNorm,1,length(checkPointNames));
figure;
bar(clusterPosRatio);
axis tight;
yl= ylim();
ylNew = [1 yl(2)];
set(gca,'yLim',ylNew,'xTick',[1:length(categoriesImmune)],'xTickLabel',categoriesImmune,'XTickLabelRotation',30);

figure;
for i=1:length(checkPointNames)
    subplot(4,1,i);
    hold on;
    for j=1:length(categoriesImmune)
        bar(j,clusterPosRatio(j,i),'FaceColor',cmap01(j,:));
    end
    axis tight;
    yl= ylim();
    ylNew = [1 yl(2)];
    set(gca,'yLim',ylNew,'xTick',[1:length(categoriesImmune)],'xTickLabel',[]);
    if i==4
        set(gca,'yLim',ylNew,'xTick',[1:length(categoriesImmune)],'xTickLabel',categoriesImmune,'XTickLabelRotation',30);
    end
    ylabel(checkPointNames(i));
    box off;
end

% plot composition per patient for each checkpoint
figure;
for i=1:length(checkPointNames)
    subplot(2,2,i);
    bar(cpPosImmuneClustersPerPatient(:,:,i),'stacked');
    colormap(cmap01);
    title (checkPointNames(i));
    axis tight;
    xlabel('Patients');
    ylabel('Number of cells');
end

% Plot number of positive cells per patient in a heatmap
pInds=[1:29,31,33:41];
checkPointNamesCo ={'LAG3','PD-1','PD-L1','IDO'};
cpTotalPosPerPatient = pData(pInds,[68,65,41,62,50]);
[sCpTotalPosPerPatient,sInds] = sortrows(cpTotalPosPerPatient);
figure;
imagesc(sCpTotalPosPerPatient);
set(gca,'xTick',1:length(checkPointNamesCo)+1,'xTickLabel',[checkPointNamesCo , 'Tregs'],'xTickLabelRotation',35, ...
    'yTick',[1:length(pInds)],'yTickLabel',pInds(sInds));
xlabel('Number of positive cells');

% Turn total counts to binary
% Threshold at 50 cells
cpTotalPosPerPatientBin = zeros(size(cpTotalPosPerPatient));
cpTotalPosPerPatientBin(cpTotalPosPerPatient>50) = 1;
[sCpTotalPosPerPatientBin,sInds] = sortrows(cpTotalPosPerPatientBin);
figure;
subplot(1,2,1);
imagesc(sCpTotalPosPerPatientBin);
set(gca,'xTick',1:length(checkPointNamesCo)+1,'xTickLabel',[checkPointNamesCo , 'Tregs'],'xTickLabelRotation',35, ...
    'yTick',[1:length(pInds)],'yTickLabel',pInds(sInds));
ylabel('Patients');
subplot(1,2,2);
totalVec = pData(pInds,38);
barh(flipud(totalVec(sInds)),'FaceColor',[0.5 0.5 0.5]);
axis tight;
set(gca,'yTick',[]);
xlabel('Number of immune cells');

% calculate conditional probability and plot
% for each checkpoint marker (e.g. pd1, calculate the probability of being
% positive for it and the probability of being positive for it given other
% markers
condProbMat = zeros(length(checkPointNamesCo)+1,length(checkPointNamesCo)+1);
for i=1:length(checkPointNamesCo)+1
    baseProb = sum(cpTotalPosPerPatientBin(:,i)) / length(cpTotalPosPerPatientBin(:,i));
    for j=1:length(checkPointNamesCo)+1
        if (i==j)
            condProbMat(i,j) = baseProb;
        else
            condProbMat(i,j) = sum(cpTotalPosPerPatientBin(:,i) & cpTotalPosPerPatientBin(:,j)) / sum(cpTotalPosPerPatientBin(:,i));
        end
    end
end

figure;
imagesc(condProbMat);
set(gca,'xTick',1:length(checkPointNamesCo)+1,'xTickLabel',[checkPointNamesCo , 'Tregs'],'xTickLabelRotation',35, ...
    'yTick',1:length(checkPointNamesCo)+1,'yTickLabel',[checkPointNamesCo , 'Tregs']);

% get the number of cells double positive for each one of the markers
markerRawInds = [16,21,23,34,15];
dpNum = zeros(length(markerRawInds),length(markerRawInds));
dataAllNo30 = dataAll;
dataAllNo30(dataAllNo30(:,1) == 30,:) = [];
for i = 1:length(markerRawInds)
    marker1 = markerRawInds(i);
    for j = 1:length(markerRawInds)
        marker2 = markerRawInds(j);
        dpNum(i,j) = sum((dataAllNo30(:,marker1) > t) & (dataAllNo30(:,marker2) > t) & (dataAllNo30(:,55) == 2));
        if (i==5) & (j~=5)
            dpNum(i,j) = sum((dataAllNo30(:,57) == 1) & (dataAllNo30(:,marker2) > t) );
        elseif (i==5) & (j==5)
            dpNum(i,j) = sum((dataAllNo30(:,57) == 1));
        elseif (i~=5) & (j==5)
            dpNum(i,j) = sum((dataAllNo30(:,57) == 1) & (dataAllNo30(:,marker1) > t) );
        end
    end
end

% Turn into proportions - divide by total number of immune cells
immuneNum = sum(dataAllNo30(:,55) == 2);
dpFrac = dpNum./immuneNum;
% look at conditional probabilities
dpCondProb = zeros(length(checkPointNamesCo)+1,length(checkPointNamesCo)+1);
for i=1:length(checkPointNamesCo)+1
    marker1 = markerRawInds(i);
    for j=1:length(checkPointNamesCo)+1
        if (i==j)
            dpCondProb(i,j) = dpFrac(i,j);
        else
            dpCondProb(i,j) = dpNum(i,j) / dpNum(i,i);
        end
    end
end

figure;
imagesc(dpCondProb);
set(gca,'xTick',1:length(checkPointNamesCo)+1,'xTickLabel',[checkPointNamesCo , 'Tregs'],'xTickLabelRotation',35, ...
    'yTick',1:length(checkPointNamesCo)+1,'yTickLabel',[checkPointNamesCo , 'Tregs']);

% draw a venn diagram of pd-1 and lag3
figure;
h=venn([dpNum(2,2),dpNum(1,1)],dpNum(1,2));
set(h(1), 'FaceColor','blue'); 
% draw a venn diagram of pd-l1 and ido
figure;
h=venn([dpNum(3,3),dpNum(4,4)],dpNum(3,4));
set(h(2), 'FaceColor',[1 0.5 0]); 
% draw a venn diagram of pd-1 and pd-l1
figure;
h=venn([dpNum(2,2),dpNum(3,3)],dpNum(2,3));
set(h(1), 'FaceColor','blue');
set(h(2), 'FaceColor','red');
% draw a venn diagram of pd-1 and ido
figure;
h=venn([dpNum(2,2),dpNum(4,4)],dpNum(2,4));
set(h(1), 'FaceColor','blue');
set(h(2), 'FaceColor',[1 0.5 0]);



% draw a venn diagram of pd-1 and lag3 and ido
% find intersection for pd1,lag3 and ido
marker1 = markerRawInds(1);
marker2 = markerRawInds(2);
marker4 = markerRawInds(4);
pd1_lag_ido = sum((dataAllNo30(:,marker1) > t) & (dataAllNo30(:,marker2) > t) & (dataAllNo30(:,marker4) > t) ...
    & (dataAllNo30(:,55) == 2));
figure;
venn([dpNum(1,1),dpNum(2,2),dpNum(4,4)],[dpNum(1,2),dpNum(1,4),dpNum(2,4),pd1_lag_ido]);

%% plot fraction of positive cells per patient in a heatmap
totalImmunePerPatientRep = repmat(pData(pInds,38),1,length(checkPointNames)+1);
totalImmunePerPatient(totalImmunePerPatientRep<50) = 0;
cpFractionPosPerPatient = cpTotalPosPerPatient ./ totalImmunePerPatientRep;
cpFractionPosPerPatient(isnan(cpFractionPosPerPatient)) = 0;
[sCpFractionPosPerPatient,sInds] = sortrows(cpFractionPosPerPatient);
figure;
imagesc(sCpFractionPosPerPatient);
set(gca,'xTick',1:length(checkPointNames)+1,'xTickLabel',[checkPointNames , 'Tregs'],'xTickLabelRotation',35, ...
    'yTick',[1:length(pInds)],'yTickLabel',pInds(sInds));
xlabel('Fraction of positive cells');

% Turn fraction to binary
% threshold at 5%
cpBinPosPerPatient = zeros(size(cpFractionPosPerPatient));
cpBinPosPerPatient(cpFractionPosPerPatient>0.1) = 1;
[sCpBinPosPerPatient,sInds] = sortrows(cpBinPosPerPatient);
figure;
subplot(1,2,1);
imagesc(sCpBinPosPerPatient);
set(gca,'xTick',1:length(checkPointNames)+1,'xTickLabel',[checkPointNames , 'Tregs'],'xTickLabelRotation',35, ...
    'yTick',[1:length(pInds)],'yTickLabel',pInds(sInds));
xlabel('Fraction of positive cells');
subplot(1,2,2)
totalVec = pData(pInds,38);
barh(flipud(totalVec(sInds)),'FaceColor',[0.5 0.5 0.5]);
axis tight;

%% Analyze lag3
% Get fraction of tumor-positive Lag3
dataTumor = dataAll((dataAll(:,55) == 5) | (dataAll(:,55) == 6),:);
tumorLag3total = sum(dataTumor(:,lag3col)>t);
tumorTotal = length(dataTumor);
percentPosTumor = tumorLag3total/tumorTotal;

% Get expression of Lag3 on different cell types
dataImmune = dataAll(dataAll(:,55) == 2,:);

% plot immune clusters for pdl-1 positive immune cells
lag3PosImmune = dataImmune(dataImmune(:,lag3col)>t,:);
lag3NegImmune = dataImmune(dataImmune(:,lag3col)<=t,:);

% for each patient get the composition of lag3-pos cells
lag3PosImmuneClustersPerPatient = zeros(points,length(categoriesImmune));
for i=[1:29,31:points]
    patientLag3PosImmune = lag3PosImmune((lag3PosImmune(:,1) == i),:);
    for j=1:length(categoriesImmune)
        lag3PosImmuneClustersPerPatient(i,j) = sum((patientLag3PosImmune(:,57) == j));
    end
end

% get clusters for all patients together
clusterPos = splitapply(@sum,ones(length(lag3PosImmune),1),lag3PosImmune(:,57));
clusterNeg = splitapply(@sum,ones(length(lag3NegImmune),1),lag3NegImmune(:,57));
clusterPosNorm = clusterPos/sum(clusterPos);
clusterNegNorm = clusterNeg/sum(clusterNeg);

% plot the fraction of LAG3-positive cells
% partition immune cells to subtypes and also plot tumor cells
figure;
PosFraction = clusterPos./(clusterPos+clusterNeg);
PosFractionAll = [PosFraction ; percentPosTumor];
categoriesAll = [categoriesImmune([1:end-1]) , 'Other Immune' , 'Tumor'];
[sVals sInds] = sort(PosFractionAll,'descend');
bar(sVals, 'FaceColor', [0.5 0.5 0.5]);
set(gca,'xTickLabel',categoriesAll(sInds),'XTickLabelRotation',20);
ylabel('Fraction of Lag3 positive cells');
axis tight;
box off;

% plot composition of LAG3-pos cells
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

% per patient
figure;
bar(lag3PosImmuneClustersPerPatient,'stacked');
colormap(cmap01);

% total
figure;
bar([clusterPosNorm, clusterNegNorm]','stacked');
colormap(cmap01);

% divide positive by negative
figure;
posDivNeg = log2(clusterPosNorm./clusterNegNorm);
[sPosDivNeg,sInds] = sort(posDivNeg,'descend');
bar(sPosDivNeg);
set(gca,'xTick',[1:length(sPosDivNeg)],'xTickLabel',categoriesImmune(sInds));
ylabel('log2[Lag3-pos / Lag3-neg]');

%% plot a bar graph of number of lag3-pos cells, pd1-cd4 and pd1-cd8
[sLag3 sInds] = sort(pData(:,5),'descend');
figure;
bar(pData(sInds,[5,10,43,44,50]));
set(gca,'xTick',[1:length(sInds)],'xTickLabel',sInds);

%% look at coocurrence. For each patient determine whether it has one of the following:
% 1. Lag3
% 2. pd1-cd4
% 3. pd1-cd8
% 4. Tregs
tBin=50;
pDataBinHeaders = {'lag3'  ; 'pd1-pos CD8' ; 'Tregs' ; 'pd1-pos CD4'};
pDataBin = zeros(points,length(pDataBinHeaders));
tmp = pData(:,[5,44,50,43]);
pDataBin(tmp>=tBin) = 1;
%[sLag3 sInds] = sort(sum(pDataBin,2),'descend');
[sPDataBin,sInds] = sortrows(pDataBin);
figure;
subplot(1,2,1);
%imagesc(pDataBin(sInds,:));
imagesc(sPDataBin);
set(gca,'xTick',[1:length(pDataBinHeaders)],'xTickLabel',pDataBinHeaders,'XTickLabelRotation',20);
set(gca,'yTick',[1:points],'yTickLabel',sInds);
subplot(1,2,2);
sIndsRev=flipud(sInds);
barh(pData(sIndsRev,38),'FaceColor',[.5 .5 .5]);
set(gca,'yTick',[]);
axis tight;

clustergram(pDataBin([1:29,31:points],:),'RowLabels',[1:29,31:points],'ColumnLabels',pDataBinHeaders,'Symmetric',false, ...
    'RowPDist','euclidean','ColumnPDist','euclidean','DisplayRatio',0.1,'ColumnLabelsRotate',20);

% look at percentage out of total:
% For PD1, look at percentange out of CD4 or CD8
% For Lag3-CD8, look at percentage out of CD8
t=20;
pDataPercentageHeaders = {'lag3-pos CD8' ; 'pd1-pos CD8' ; 'pd1-pos CD4'};
pDataPercentage = zeros(points,length(pDataPercentageHeaders));
a = pData(:,49);
b = pData(:,53);
a(a<t) = 0;
b(a<t) = 0;
pDataPercentage(:,1) = a./b;
a = pData(:,44);
b = pData(:,53);
a(a<t) = 0;
b(a<t) = 0;
pDataPercentage(:,2) = a./b;
a = pData(:,43);
b = pData(:,52);
a(a<t) = 0;
b(a<t) = 0;
pDataPercentage(:,3) = a./b;
pDataPercentage(isnan(pDataPercentage)) = 0;
[sLag3, sInds] = sort(pDataPercentage(:,1),'descend');
figure;
bar(pDataPercentage(sInds,:));
set(gca,'xTick',[1:points],'xTickLabel',sInds);
legend(pDataPercentageHeaders);

%% Examine coocurrence of lag3-pos and pd1-pos cd8 cells
pDataCooccurenceHeaders = {'% lag3-pos CD8' ; '% pd1-pos CD8' ; 'expected cooccurence' ; '% lag3-pos pd1-pos cd8'};
pDataCooccurence = pDataPercentage(:,[1:2]);
pDataCooccurence(:,3) = pDataCooccurence(:,1) .* pDataCooccurence(:,2);
a = pData(:,54);
b = pData(:,53);
a(a<t) = 0;
b(a<t) = 0;
pDataCooccurence(:,4) = a./b;

[sPDataCooccurence, sInds] = sortrows(pDataCooccurence,'descend');
figure;
bar(sPDataCooccurence);
set(gca,'xTick',[1:points],'xTickLabel',sInds);
legend(pDataCooccurenceHeaders);
figure;
imagesc(sPDataCooccurence);
set(gca,'yTick',[1:points],'yTickLabel',sInds);
set(gca, 'xTick',[1:length(pDataCooccurenceHeaders)],'xTickLabel',pDataCooccurenceHeaders);

%% Examine coocurrence of lag3-pos and pd1-pos cd4 cells
pDataCooccurenceHeaders = {'% lag3-pos CD4' ; '% pd1-pos CD4' ; 'expected cooccurence' ; '% lag3-pos pd1-pos cd4'};
a = pData(:,48);
b = pData(:,52);
a(a<t) = 0;
b(a<t) = 0;
pDataCooccurence(:,1) = a./b;
pDataCooccurence(:,2) = pDataPercentage(:,3);
pDataCooccurence(:,3) = pDataCooccurence(:,1) .* pDataCooccurence(:,2);
a = pData(:,55);
b = pData(:,52);
a(a<t) = 0;
b(a<t) = 0;
pDataCooccurence(:,4) = a./b;
pDataCooccurence(isnan(pDataCooccurence)) = 0;
[sPDataCooccurence, sInds] = sortrows(pDataCooccurence,'descend');
figure;
bar(sPDataCooccurence);
set(gca,'xTick',[1:points],'xTickLabel',sInds);
legend(pDataCooccurenceHeaders);
figure;
imagesc(sPDataCooccurence);
set(gca,'yTick',[1:points],'yTickLabel',sInds);
set(gca, 'xTick',[1:length(pDataCooccurenceHeaders)],'xTickLabel',pDataCooccurenceHeaders);

%% Examine coocurrence of lag3-pos and pd1-pos T cells
pDataCooccurenceHeaders = {'% lag3-pos T cells' ; '% pd1-pos Tcells' ; 'expected cooccurence' ; '% lag3-pos pd1-pos Tcells'};
a = pData(:,57);
b = pData(:,59);
a(a<t) = 0;
b(a<t) = 0;
pDataCooccurence(:,1) = a./b;
a = pData(:,56);
b = pData(:,59);
a(a<t) = 0;
b(a<t) = 0;
pDataCooccurence(:,2) = a./b;
pDataCooccurence(:,3) = pDataCooccurence(:,1) .* pDataCooccurence(:,2);
a = pData(:,58);
b = pData(:,59);
a(a<t) = 0;
b(a<t) = 0;
pDataCooccurence(:,4) = a./b;
pDataCooccurence(isnan(pDataCooccurence)) = 0;
[sPDataCooccurence, sInds] = sortrows(pDataCooccurence,'descend');
figure;
bar(sPDataCooccurence);
set(gca,'xTick',[1:points],'xTickLabel',sInds);
legend(pDataCooccurenceHeaders);
figure;
imagesc(sPDataCooccurence);
set(gca,'yTick',[1:points],'yTickLabel',sInds);
set(gca, 'xTick',[1:length(pDataCooccurenceHeaders)],'xTickLabel',pDataCooccurenceHeaders);

%%
% plot a bar graph of number of Tregs
figure;
[tregNumSort, sInds] = sort(pData(:,51)./pData(:,50));
bar(tregNumSort);
set(gca,'xTick',1:length(sInds),'xTickLabel',sInds);

% %% cell scatter
% figure;
% ind1 = 16;
% ind2 = 21;
% scatter(dataAll(:,ind1), dataAll(:,ind2));
% xlabel(dataHeaders{ind1});
% ylabel(dataHeaders{ind2});

%% test graphs:

% dot plots of patients - lag3 vs. pd1
ind1 = 5;
ind2 = 10;
figure;
subplot(3,1,1);
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
for i = 1:points
    text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_lag3_pd1,PVAL_lag3_pd1] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1));

% dot plots of patients -  tregs vs. pd1
ind1 = 50;
ind2 = 10;
subplot(3,1,2);
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
for i = 1:points
    text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_treg_pd1,PVAL_treg_pd1] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1))

% dot plots of patients -  tregs vs. lag3
ind1 = 50;
ind2 = 5;
subplot(3,1,3);
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
for i = 1:points
    text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_treg_lag3,PVAL_treg_lag3] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1))

%% good graphs:
% 1. # of pdl-1 positive immune cells vs. pd1. positive corr.
% 2. # of pdl-1 positive tumor cells vs. ratio of CD8 and CD4 pd1 

% dot plots of patients - pd-l1 vs. pd1
ind1 = 12;
ind2 = 10;
figure;
subplot(3,1,1);
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
% for i = 1:points
%     text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
% end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_pdl1_pd1,PVAL_pdl1_pd1] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1));


% dot plots of patients -  immune pd-l1 vs. pd1
ind1 = 41;
ind2 = 10;
subplot(3,1,2);
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
% for i = 1:points
%     text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
% end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_pdl1Immune_pd1,PVAL_pdl1Immune_pd1] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1))

% dot plots of patients - pdl1 tumor vs. pd1 cd4-cd8
ind1 = 40;
ind2 = 45;
% plot log
subplot(3,1,3);
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2));
scatter(X,Y,36,[0 0 0],'filled');
% for i = 1:points
%     text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
% end
% p = polyfit(X,Y,2);
% f = polyval(p,X);
% [sX sInds] = sort(X);
% plot(sX,f(sInds),'k');
[fitpoly2 gof]=fit(X,Y,'poly2');
gof
Ypred = fitpoly2(X);
[sX sInds] = sort(X);
plot(sX,Ypred(sInds),'Color',[0.7 0.7 0.7]);

xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_pdl1Tumor_pd1Cd4Cd8,PVAL_pdl1Tumor_pd1Cd4Cd8] = corr(X,Y)

xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));

% dot plots of patients -  tumor pd-l1 vs. pd1
ind1 = 40;
ind2 = 10;
figure
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
% for i = 1:points
%     text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
% end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_pdl1Tumor_pd1,PVAL_pdl1Tumor_pd1] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1))

% dot plots of patients -  immune cells vs. pd1
ind1 = 38;
ind2 = 10;
figure
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2)+1);
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
% for i = 1:points
%     text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
% end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_immuneNum_pd1,PVAL_immuneNum_pd1] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1))

% dot plots of patients -  immune cells vs. cd4-8 ratio
ind1 = 38;
ind2 = 45;
figure
hold on;
X = log2(pData(:,ind1)+1);
Y = log2(pData(:,ind2));
scatter(X,Y,36,[0.2 0.2 0.2],'filled');
% for i = 1:points
%     text(log2(pData(i,ind1)+1)+0.1,log2(pData(i,ind2)+1),num2str(i));
% end
p = polyfit(X,Y,1);
f = polyval(p,X);
plot(X,f,'Color',[0.7 0.7 0.7]);
xlabel(pHeaders(ind1));
ylabel(pHeaders(ind2));
axis tight;
[RHO_immuneNum_Cd4Cd8,PVAL_immuneNum_Cd4Cd8] = corr(log2(pData(:,ind1)+1),log2(pData(:,ind2)+1))

% bar plots of CD4/CD8 pd1 ratio

% dot plot of cells
ind1 = 16;
ind2 = 21;
ind3 = 32;
figure;
scatter3(dataAll(:,ind1),dataAll(:,ind2),dataAll(:,ind3));
xlabel(dataHeaders(ind1));
ylabel(dataHeaders(ind2));
zlabel(dataHeaders(ind3));

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
tumorLag3total = sum(pData(:,4));
immunePD1total = sum(pData(:,5));
immuneTotal = sum(pData(:,1));
tumorTotal = sum((dataAll(:,55) == 5) | (dataAll(:,55) == 6));
percentPosImmune = immunePD1total/immuneTotal;
percentPosTumor = tumorLag3total/tumorTotal;
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
    patientLag3PosImmune = pdl1PosImmune((pdl1PosImmune(:,1) == i),:);
    for j=1:length(categoriesImmune)
        pdl1PosImmuneClustersPerPatient(i,j) = sum((patientLag3PosImmune(:,57) == j));
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
figure;
PosFraction = clusterPos./(clusterPos+clusterNeg);
[sVals sInds] = sort(PosFraction,'descend');
bar(sVals, 'FaceColor', [0.5 0.5 0.5]);
set(gca,'xTickLabel',categoriesImmune(sInds),'XTickLabelRotation',20);
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
