% Mibi normalize batches of panCK
% We had two staining sessions. One probably had more antibody for panCK
% because all cores in that batch have higher panCK. Need to normalize.
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
vec=1:coreNum;
channelInd = 42;
imsize = 2048;
batch1Inds = [1:29];
batch2Inds = [30:41];
batch1SampleNum = length(batch1Inds);
batch2SampleNum = length(batch2Inds);

% %load data
dataPerChannel = zeros(imsize,imsize,coreNum);
dataReshape = zeros(imsize*imsize,coreNum);
for i=vec
    load([corePath{i},'dataDeNoise.mat']);
    currData=countsNoBg(:,:,channelInd);
    dataPerChannel(:,:,i) = currData;
    dataReshape(:,i) = currData(:);
end
disp('finished loading');

% create a vector of batch1 and of batch2
batch1 = dataReshape(:,batch1Inds);
batch1 = reshape(batch1,imsize*imsize*batch1SampleNum,1);

batch2 = dataReshape(:,batch2Inds);
batch2 = reshape(batch2,imsize*imsize*batch2SampleNum,1);

% subsample batch1 to the length of batch 2
batch1subSample = datasample(batch1,length(batch2),'Replace',false);

% plot histograms
figure;
hold on;
histogram(batch1subSample,[1:1:30]);
histogram(batch2,[1:1:30]);

% normalize batch2 to batch 1 by quantile normalization
% we sort both array to create a mapping function. Each value in batch2
% will be mapped to the corresponding value in batch 1. Note that some
% values in batch2 will be mapped to 2 values in batch 1.
batch1S = sort(batch1subSample);
batch2S = sort(batch2);
figure;
hold on;
plot(batch1S(1:100:end),batch2S(1:100:end));
xlabel('batch 1');
ylabel('batch 2');
xl=xlim();
yl=ylim();
maxVal = max([xl,yl]);
plot([0,maxVal],[0,maxVal]);
% remove overlapping values
[batch2SU,ia,ic] = unique(batch2S,'legacy');
batch1SU = batch1S(ia);
figure; 
hold on;
plot(batch1SU,batch2SU);
xlabel('batch 1');
ylabel('batch 2');
xl=xlim();
yl=ylim();
maxVal = max([xl,yl]);
plot([0,maxVal],[0,maxVal]);

% 
% for each sample in batch2, transform the data and write to file
for i=1:length(batch2Inds)
    disp(i);
    sampleID = batch2Inds(i);
    currData = batch2([(i-1)*(imsize*imsize)+1:i*(imsize*imsize)]);
    % transform
    [tf, Inds] = ismember(currData,batch2SU);
    currDataTrans = batch1SU(Inds);
    currDataTransMat = reshape(currDataTrans,imsize,imsize);
    % save
    load([corePath{sampleID},'dataDeNoise.mat']);
    countsNoBg(:,:,channelInd) = currDataTransMat;
    save ([corePath{sampleID},'dataDeNoise.mat'],'totalIonFiltSum','countsAllSFiltCRSum','countsNoNoise','massDS','IntNormD','pointNumber','countsNoBg');
end




% % plot intensity distributions for both batches
% hedges = [1:1:50];
% hline=zeros(coreNum,length(hedges)-1);
% for j = vec
%     data = p{j}.countsNoBg(:,:,plotChannelInd);
%     h=histogram(data(:),hedges);
%     %h=histogram(data,hedges,'Normalization','cdf');
%     hline(j,:)=h.Values;
% end
% clear('h','data');
% a = 1:coreNum ;
% labels = strread(num2str(a),'%s');
% figure;
% plot(hedges([1:end-1]),hline,'-*');
% legend(labels);
