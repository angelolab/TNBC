% check background removal

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
    '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/Point19/' };

coreNum= length(corePath);

%load data
p=cell(coreNum,1);
for i=[1:15]
    i
    p{i}=load([corePath{i},'dataDeNoise.mat']);
    p{i}.IntNormD50={};
    p{i}.countsAllSFiltCRSumGaus=[];
    p{i}.countsNoNoise=[];
    p{i}.totalIonFiltSum=[];
    if i>1
        p{i}.massDS=[];
    end
end

% %for each print the NN dist
% plotChannelInd=10;
% for j = [1,13,21:28]
%     p{j}.IntNormDPre{plotChannelInd}=MibiGetIntNormDist(p{j}.countsAllSFiltCRSum(:,:,plotChannelInd),p{j}.countsAllSFiltCRSum(:,:,plotChannelInd),50,2,50);
% end
for j = [1,13,21:28]
    mask(:,:,j) = MibiGetMask(p{j}.countsAllSFiltCRSum(:,:,10),cap);
end
    
f=figure;
hedges = [0:0.25:30];
hline=zeros(coreNum,length(hedges)-1);
for j = [1,13,21:28]
    data = p{j}.IntNormDPre{plotChannelInd};
    h=histogram(data,hedges,'Normalization','probability');
    %h=histogram(data,hedges,'Normalization','cdf');
    hline(j,:)=h.Values;
end
clear('h','data');
a = 1:coreNum ;
labels = strread(num2str(a),'%s');
plot(hedges([1:end-1]),hline);
legend(labels);

countsNoNoise = zeros(size(p{j}.countsNoBg,1),size(p{j}.countsNoBg,2),coreNum);
for i=[1,13,21:28]
    countsNoNoise(:,:,i) = MibiFilterImageByNNThreshold(p{i}.countsAllSFiltCRSum(:,:,10),p{i}.IntNormDPre{10},5);
end

% for each print the background channel before and after
% noise removal
capImage=5;
for i=[1,13,21:28]
    for j=[7,12,10]
        % plot before
        data = p{i}.countsAllSFiltCRSum(:,:,j);
        data(data>capImage)=capImage;
        figure;
        imagesc(data);
        title(['Before - point',num2str(i),' - ', p{1}.massDS.Label(j)]);
        plotbrowser on;
        % plot after using method 1
        data = p{i}.countsNoBg(:,:,j);
        data(data>capImage)=capImage;
        figure;
        imagesc(data);
        title(['After 1 - point',num2str(i),' - ',p{1}.massDS.Label(j)]);
        plotbrowser on;
    end
%     for j=10
%         % plot after
%         data = countsNoNoise(:,:,i);
%         data(data>capImage)=capImage;
%         figure;
%         imagesc(data);
%         title(['After 2 - point',num2str(i),' - ',p{1}.massDS.Label(j)]);
%         plotbrowser on;
%     end
end


% % for each print the dna,background and foxp3 channels before and after
% % noise removal
% capImage=20;
% for i=[1,13,21:28]
%     for j=10
%     %for i=21
%     %for j=[10,7,12]
%         % plot before
%         data = p{i}.countsAllSFiltCRSum(:,:,j);
%         data(data>capImage)=capImage;
%         figure;
%         imagesc(data);
%         title(['Before - point',num2str(i),' - ', p{1}.massDS.Label(j)]);
%         plotbrowser on;
%         % plot after
%         data = p{i}.countsNoBg(:,:,j);
%         data(data>capImage)=capImage;
%         figure;
%         imagesc(data);
%         title(['After - point',num2str(i),' - ',p{1}.massDS.Label(j)]);
%         plotbrowser on;
%     end
% end