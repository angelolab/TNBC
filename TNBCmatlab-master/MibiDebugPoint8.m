% % Debug patient8
% points=41;
% massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
% path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
% resultsDir = [path,'/dataPerCell'];
% 
% % add two cols at the beginning for patient index and label vec
% dataCellsAll = zeros(1,51);
% dataCellsTransAll = zeros(1,51);
% dataScaleSizeCellsAll = zeros(1,51);
% dataScaleSizeCellsTransAll = zeros(1,51);
% labelVecAll = zeros(1,1);
% 
% %% load data
% insertInd = 1;
% for i=1:points
%     load([path,'/Point',num2str(i),'/cellData.mat']);
%     currLength = length(dataCells);
%     % load data
%     dataCellsAll([insertInd:insertInd+currLength-1],[3:end]) = dataCells;
%     dataCellsTransAll([insertInd:insertInd+currLength-1],[3:end]) = dataCellsTrans;
%     dataScaleSizeCellsAll([insertInd:insertInd+currLength-1],[3:end]) = dataScaleSizeCells;
%     dataScaleSizeCellsTransAll([insertInd:insertInd+currLength-1],[3:end]) = dataScaleSizeCellsTrans;
%     % add label
%     dataCellsAll([insertInd:insertInd+currLength-1],2) = labelVec;
%     dataCellsTransAll([insertInd:insertInd+currLength-1],2) = labelVec;
%     dataScaleSizeCellsAll([insertInd:insertInd+currLength-1],2) = labelVec;
%     dataScaleSizeCellsTransAll([insertInd:insertInd+currLength-1],2) = labelVec;
%     % add patient index
%     dataCellsAll([insertInd:insertInd+currLength-1],1) = i;
%     dataCellsTransAll([insertInd:insertInd+currLength-1],1) = i;
%     dataScaleSizeCellsAll([insertInd:insertInd+currLength-1],1) = i;
%     dataScaleSizeCellsTransAll([insertInd:insertInd+currLength-1],1) = i;
%     labelVecAll([insertInd:insertInd+currLength-1]') = labelVec;
%     insertInd = insertInd+currLength;
% end

% %% standardize
% % standardize either by quantile normalization or std normalization
% dataScaleSizeCellsTransAllStd = zscore(dataScaleSizeCellsTransAll);
% dataCellsTransAllStd = zscore(dataCellsTransAll);
% dataCellsAllStd = zscore(dataCellsAll);
% dataScaleSizeCellsAllStd = zscore(dataScaleSizeCellsAll);
% 
% dataScaleSizeCellsTransAllStd(:,[1:2]) = dataScaleSizeCellsTransAll(:,[1:2]);
% dataCellsTransAllStd(:,[1:2]) = dataCellsTransAll(:,[1:2]);
% dataCellsAllStd(:,[1:2]) = dataCellsAll(:,[1:2]);
% dataScaleSizeCellsAllStd(:,[1:2]) = dataScaleSizeCellsAll(:,[1:2]);

% for all patients plot a histogram of their channel X values in table Y
channelInd = 39;
hedges = [0:10:200];
hline=zeros(points,length(hedges)-1);
for i=1:points
    pInds = (dataCellsAll(:,1)==i);
    currData = dataCellsAll(pInds,channelInd);
    h=histogram(currData,hedges,'Normalization','probability');
    hline(i,:)=h.Values;
end

a = 1:points ;
labels = strread(num2str(a),'%s');
figure;
plot(hedges([1:end-1]),hline);
legend(labels);
