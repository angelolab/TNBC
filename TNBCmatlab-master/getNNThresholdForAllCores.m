% get nn-filter threshold for all cores

corePath = {'/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA460-S1R8C2/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA460-S1R4C5/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA460-S1R6C9/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA460-S2R3C3/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA460-S2R9C5/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_MultipleCores/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_MultipleCores2/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_MultipleCores3/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores/Point2/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point2/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point3/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point4/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point5/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point6/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point7/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point8/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point9/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point10/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point11/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point12/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point13/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point14/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point15/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point16/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point17/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point18/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point19/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining/Point2/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point1/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point2/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point3/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point4/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point5/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point6/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point7/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point8/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point9/', ...
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-460_secondstaining2/Point10/'};

coreNum= length(corePath);
vec=[1:2,4,9:10,16,30:35,39:40];
% % load data
% p=cell(coreNum,1);
% for i=vec
%     i
%     p{i}=load([corePath{i},'dataDeNoise.mat']);
%     p{i}.countsAllSFiltCRSum=[];
%     p{i}.countsAllSFiltCRSumGaus=[];
%     p{i}.countsNoNoise=[];
%     p{i}.totalIonFiltSum=[];
%     if i>1
%         p{i}.massDS=[];
%     end
% end

disp('finished loading');

% plot a channel from all cores
plotChannelInd=39;
capImage = 5;

% get the NN values for the channel for all cores
% for i=vec
%     p{i}.IntNormD{plotChannelInd}=MibiGetIntNormDist(p{i}.countsNoBg(:,:,plotChannelInd),p{i}.countsNoBg(:,:,plotChannelInd),25,2,25);
% end

% chanelInAllPoints = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);
% chanelInAllPointsCapped = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);

% % plot individual plots
% for j = 1:coreNum
%     chanelInAllPoints(:,:,1,j) = p{j}.countsNoBg(:,:,plotChannelInd);
%     data = chanelInAllPoints(:,:,1,j);
%     data(data>capImage) = capImage;
%     chanelInAllPointsCapped(:,:,1,j) = data;
%     figure;
%     imagesc(data);
%     title(['Core ',num2str(j)]);
%     plotbrowser on;
% end
% plot together
%f1=figure;
%imdisp(chanelInAllPointsCapped,'Size',[4,7],'Border',[0.005,0.005]);

% % plot the NN histograms for all points
% f=figure;
% [ha, pos] = tight_subplot(4,ceil(coreNum/4),[.05 .05],[.05 .05],[.05 .05]); 
% for j = 1:coreNum
%     axes(ha(j));
%     data = p{j}.IntNormD{plotChannelInd};
%     histogram(data,100,'Normalization','probability');
%     set(ha,'xlim',[0 30]);
% end

% plot the NN histograms for all points in a single plot
f=figure;
hedges = [0:0.25:30];
hline=zeros(coreNum,length(hedges)-1);
for j = vec
    data = p{j}.IntNormD{plotChannelInd};
    h=histogram(data,hedges,'Normalization','probability');
    %h=histogram(data,hedges,'Normalization','cdf');
    hline(j,:)=h.Values;
end
clear('h','data');
a = 1:coreNum ;
labels = strread(num2str(a),'%s');
plot(hedges([1:end-1]),hline);
legend(labels);

% %% scale the histograms so that they will have the same 1 and 95 percentile
% 
% % scale both sides
% prc95 = zeros(coreNum,1);
% prc1 = zeros(coreNum,1);
% setVal95 = 0.95;
% setVal1 = 0.01;
% 
% for j=1:coreNum
%     data = p{j}.IntNormD{plotChannelInd};
%     % scale data between 0 and 1
%     maxVal = max(data);
%     minVal = min(data);
%     data01 = (data-minVal)./maxVal;
%     prc95(j) = prctile(data01,95);
%     prc1(j) = prctile(data01,1);
%     dataScaled = imadjust(data01,[prc1(j); prc95(j)],[setVal1; setVal95]);
%     %dataScaled = data*(setVal95/prc95(j));
%     p{j}.IntNormDscaled{plotChannelInd}=dataScaled;
% end

% % scale high side
% % setVal95 = 20;
% % 
% % for j=1:coreNum
% %     data = p{j}.IntNormD{plotChannelInd};
% %     prc95(j) = prctile(data,95);
% %     dataScaled = data*(setVal95/prc95(j));
% %     p{j}.IntNormDscaled{plotChannelInd}=dataScaled;
% % end

% % plot the scaled NN histograms for all points in a single plot
% f=figure;
% hedges = [0:0.01:1];
% hlineScaled=zeros(coreNum,length(hedges)-1);
% for j = 1:coreNum
%     data = p{j}.IntNormDscaled{plotChannelInd};
%     h=histogram(data,hedges,'Normalization','probability');
%     %h=histogram(data,hedges,'Normalization','cdf');
%     hlineScaled(j,:)=h.Values;
% end
% clear('h','data');
% a = 1:coreNum ;
% labels = strread(num2str(a),'%s');
% plot(hedges([1:end-1]),hlineScaled);
% legend(labels);

% % plot the NN histograms distance 50 for all points in a single plot
% f=figure;
% hedges = [0:0.25:30];
% hline=zeros(coreNum,length(hedges)-1);
% for j = 1:coreNum
%     data = p{j}.IntNormD50{plotChannelInd};
%     h=histogram(data,hedges,'Normalization','probability');
%     hline(j,:)=h.Values;
% end
% a = 1:coreNum ;
% labels = strread(num2str(a),'%s');
% clear('h','data');
% plot(hedges([1:end-1]),hline);
% legend(labels);

% remove noise and replot
countsNoNoise = zeros(size(p{j}.countsNoBg,1),size(p{j}.countsNoBg,2),coreNum);
for j = vec
    countsNoNoise(:,:,j) = MibiFilterImageByNNThreshold(p{j}.countsNoBg(:,:,plotChannelInd),p{j}.IntNormD{plotChannelInd},3);
end

% % replot
% f3=figure;
% chanelNoNoiseInAllPoints = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);
% chanelNoNoiseInAllPointsCapped = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);
% 
% for j = 1:coreNum
%     chanelNoNoiseInAllPoints(:,:,1,j) = mat2gray(countsNoNoise(:,:,j));
%     data = chanelNoNoiseInAllPoints(:,:,1,j);
%     data(data>capImage) = capImage;
%     chanelNoNoiseInAllPointsCapped(:,:,1,j) = data;
% end
% imdisp(chanelNoNoiseInAllPointsCapped,'Size',[4,7],'Border',[0.005,0.005]);
% clear('data');

% gaus and replot
chanelNoNoiseInAllPointsGaus = zeros(size(p{1}.countsNoBg(:,:,plotChannelInd),1),size(p{1}.countsNoBg(:,:,plotChannelInd),2),1,coreNum);

for j = vec
    % plot before
    data = p{j}.countsNoBg(:,:,plotChannelInd);
    data(data>capImage) = capImage;
    figure;
    imagesc(data);
    title(['Core ',num2str(j)]);
    plotbrowser on;
    % plot after
    data = countsNoNoise(:,:,j);
    data(data>capImage) = capImage;
    dataGaus = imgaussfilt(data,1);
    chanelNoNoiseInAllPointsGaus(:,:,1,j) = mat2gray(dataGaus);
    figure;
    imagesc(dataGaus);
    title(['Core denoise ',num2str(j)]);
    plotbrowser on;
end
%f4=figure;
%imdisp(chanelNoNoiseInAllPointsGaus,'Size',[4,7],'Border',[0.005,0.005]);
clear('data','dataGaus');