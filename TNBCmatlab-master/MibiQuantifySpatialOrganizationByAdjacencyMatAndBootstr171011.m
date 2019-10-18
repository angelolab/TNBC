% % MibiQuantifySapatialOrganizationByAdjacencyMatAndBootstrap
% Get enrichment of specific cell types to sit together or not
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
load([pathSegment,'/dataWithTumorAndImmuneGroups170911.mat']);

%% redefine cd3-positive B cells as cd4 t-cells
%figure; histogram(dataAll(dataAll(:,57)==6,33),1000);
tcd3 = 0.5;
dataAll( ((dataAll(:,57) == 6) & (dataAll(:,33) > tcd3)),57) = 2;

%% cell_types to compare:
% Go over all pairwise combinations of markers
markerInds = [11,12,15:19,21:28,30:50];
dataMarkers = dataAll(:,markerInds);
markerTitles = dataHeaders(markerInds);
markerNum = length(markerTitles);

BootstrapNum = 100;
distLim = 100;
closeNum = zeros(markerNum,markerNum);
closeNumRand = zeros(markerNum,markerNum,BootstrapNum);

t = 1;
% for i=2:41
%     disp(['Working on point:',num2str(i)]);
%     load([pathSegment,'/Point',num2str(i),'/cellDistances.mat']);
%     % get data for current patient
%     patientInds = (dataAll(:,1) == i);
%     patientData = dataAll(patientInds,:);
%     patientDataMarkers = dataMarkers(patientInds,:);
%     % go over markers
%     parfor j=1:markerNum
%         j
%         marker1PosInds = ( patientDataMarkers(:,j) > t );
%         marker1PosLabels = patientData(marker1PosInds,2);
%         marker1Num = length(marker1PosLabels);
%         for k=1:markerNum
%             marker2PosInds = ( patientDataMarkers(:,k) > t );
%             marker2PosLabels = patientData(marker2PosInds,2);
%             marker2Num = length(marker2PosLabels);
%             truncDistMat = distancesMat(marker1PosLabels,marker2PosLabels);
%             % turn to binary
%             truncDistMatBin = zeros(size(truncDistMat));
%             truncDistMatBin(truncDistMat<distLim) = 1;
%             % record interaction num
%             closeNum(j,k) = sum(sum(truncDistMatBin));
%             % randomize
%             for r=1:BootstrapNum
%                 % get random labels for marker 1
%                 marker1LabelsRand = datasample(patientData(:,2),marker1Num);
%                 % get random labels for marker 2
%                 marker2LabelsRand = datasample(patientData(:,2),marker2Num);
%                 randTruncDistMat = distancesMat(marker1LabelsRand,marker2LabelsRand);
%                 % turn to binary
%                 randTruncDistMatBin = zeros(size(randTruncDistMat));
%                 randTruncDistMatBin(randTruncDistMat<distLim) = 1;
%                 % record interaction num
%                 closeNumRand(j,k,r) = sum(sum(randTruncDistMatBin));
%             end
%         end
%     end
% 
%     % var = immuneTumorInteractionsB./(immuneInteractionsB+immuneTumorInteractionsB);
%     % realVar = immuneTumorInteractions./(immuneInteractions+immuneTumorInteractions);
%     z = zeros(markerNum);
%     muhat = zeros(markerNum);
%     sigmahat = zeros(markerNum);
% 
%     for j=1:markerNum
%         for k=1:markerNum
%             tmp= reshape(closeNumRand(j,k,:),BootstrapNum,1);
%             [muhat(j,k),sigmahat(j,k)] = normfit(tmp);
%             z(j,k) = (muhat(j,k)-closeNum(j,k))/sigmahat(j,k);
%         end
%     end
%     
%     % make matrix symmetric
%     for j=1:markerNum
%         for k=1:j
%             z(j,k) = z(k,j); 
%         end
%     end
% %     zplot = -z;
% %     zplot(isinf(zplot)) = 0;
% %     clustergram(zplot([1:15,17:end],[1:15,17:end]),'RowLabels', markerTitles([1:15,17:end]), 'ColumnLabels', markerTitles([1:15,17:end]), ...
% %         'Colormap', 'redbluecmap','DisplayRange', 50, 'DisplayRatio', 0.1);
% 
%     % save results
%     save([pathSegment,'/Point',num2str(i),'/spatialAnalysis.mat'],'closeNum','closeNumRand','z','muhat','sigmahat');
% 
% end


% visualize
i=6;
load([pathSegment,'/Point',num2str(i),'/spatialAnalysis.mat']);
zplot = -z;
zplot(isnan(zplot)) = 0;
zplot(isinf(zplot)) = 0;
clustergram(zplot([1:15,17:end],[1:15,17:end]),'RowLabels', markerTitles([1:15,17:end]), 'ColumnLabels', markerTitles([1:15,17:end]), ...
    'Colormap', 'redbluecmap','DisplayRange', 50, 'DisplayRatio', 0.1);

