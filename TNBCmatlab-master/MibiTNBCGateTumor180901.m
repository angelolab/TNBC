% MibiTNBCMarkTumorClusters170828
% for all cells add a column marking their 'primary' population
points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/';
pathCluster = [path,'ClusteringResults/'];
resultsDir = [path,'dataPerCell170809/ClusterAll170831/GateTumor170901/'];

phenoS = importdata([pathCluster,'/data_clusterTumor170831.csv']);

markersToGate = {' EGFR',' Keratin17',' Keratin6',' p53',' Pan_Keratin'};
markerCols = find(ismember(phenoS.colheaders,markersToGate));
Tgate = ([-0.3,-0.265,-0.282,-0.597,-0.535]);
% plot histograms for all markers
% for i=1: length(markersToGate)
%     figure;
%     currCol=markerCols(i);
%     histogram(phenoS.data(:,currCol),10000);
%     title(markersToGate{i});
% end

% mark cells which are positive for all markers
tumorInds = zeros(size(phenoS.data,1),1);
for i=1:length(markersToGate)
    tumorInds(phenoS.data(:,markerCols(i)) > Tgate(i)) = 1;
end

sum(tumorInds);

%% visualize tumor cells on specific patient
cyt_session_name = 'ClusteringResults/data_clusterTumor170831.mat';
sampleName = 'dataStdTumor_p26';
pointNumber=26;
gateNum=1;
phenoCol=70;

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentPerim';

colormap('parula');
load([path,'/Point',num2str(pointNumber),'/dataDeNoiseCohortNoAgg.mat']);
load([pathSegment,'/Point',num2str(pointNumber),'/segmentationParams.mat']);
load([pathSegment,'/Point',num2str(pointNumber),'/cellData.mat']);

%% get data from gate in cyt session
load(cyt_session_name);
sessionData(:,phenoCol) = tumorInds;

% get all data from sample in cyt session
[~,sampleInd] = ismember(sampleName,gates(:,1));
dataSampleInd = gates{sampleInd,2};
cellsInSample = sessionData(dataSampleInd,:);

% get points in gate
cellsInGate = cellsInSample(cellsInSample(:,phenoCol)==gateNum,:);

% plot nuclei with perim
%for i=27
for i=[7,8,9,40,37,42,46,36,21,19]
    maxv=5;
    if i==7
        maxv=25;
    end
    rgb_image = MibiGetRGBimageFromMat(countsNoNoiseNoAgg(:,:,i),maxv);
    rgb_image_perim= imoverlay(rgb_image , cellPerimNewMod , [1 .3 .3]);
    figure; imagesc(rgb_image_perim);
    title(massDS.Label{i});
    plotbrowser on;
end

%% assign populations to the dataset
stats = regionprops(newLmod,'PixelIdxList');
labelNum = length(labelIdentityNew2);
imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
counts=0;
for i=1:labelNum-1
    if labelIdentityNew2(i) == 0
        imageL(stats(i).PixelIdxList)=1;
    elseif ismember(i,cellsInGate(:,1))
        imageL(stats(i).PixelIdxList)=3;
        counts=counts+1;
    elseif ismember(i,cellsInSample(:,1))
        imageL(stats(i).PixelIdxList)=2;
    else
        imageL(stats(i).PixelIdxList)=4;
    end
end

imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));

figure;
colormap('parula');
imagesc(label2rgb(imageLReshape));
plotbrowser on;

%% save classification
dataMat = [phenoS.data(:,[2:end]) , tumorInds];
Labels2Write = [ phenoS.colheaders([2:end]) , 'tumorYN' ];

% convert values in 1st colums to point number
% cyt saves the order sorted as a string, so the numbers don't represent
% our points
pointStr = cellfun(@num2str,num2cell(1:points),'uniformoutput',0);
pointStrSort=sort(pointStr);
pointsByStrCell = cellfun(@str2num,pointStrSort,'uniformoutput',0);
pointsByStr = cell2mat(pointsByStrCell);

newPointVec = zeros(length(phenoS.data),1);
for i=1:points
    pointFromCyt = i;
    realPoint = pointsByStr(i);
    newPointVec(phenoS.data(:,1) == pointFromCyt) = realPoint;
end

%for each cell add a category label and save fcs files and data structures
%for all cells, for immune and for tumor
mkdir(resultsDir);
mkdir([resultsDir,'/dataStdRest']);
mkdir([resultsDir,'/dataStdTumor']);
TEXT.PnS = Labels2Write;
TEXT.PnN = Labels2Write;

for i=1:points
    dataPoint = dataMat((newPointVec == i),:);
    dataRest = dataPoint((dataPoint(:,end) == 0),:);
    dataTumor = dataPoint((dataPoint(:,end) == 1),:);
    writeFCS([resultsDir,'/dataStdRest/dataStdRest_p',num2str(i),'.fcs'],dataRest,TEXT); 
    writeFCS([resultsDir,'/dataStdTumor/dataStdTumor_p',num2str(i),'.fcs'],dataTumor,TEXT);
end

% generate a struct for all the data to save
phenoSOld = phenoS;
phenoS.data = [newPointVec , dataMat];
phenoS.textdata = ['GateSource' , Labels2Write];
phenoS.colheaders = ['GateSource' , Labels2Write];

save([pathCluster,'dataGateTumor170901.mat'],'phenoS');