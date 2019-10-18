corePath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/';
coreNum = 25;

% % load data
% massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
% p=cell(coreNum,1);
% for i=[1:coreNum]
%     p{i}=load([corePath,'/Point',num2str(i),'/dataDeNoiseCohort.mat']);
% end

% % parameters for filtering noise in B7H3
currChannel = 'HLA_Class_1';
gausRad=1;
gausFlag=1;
t=150;
[~,currChannelInd] = ismember(currChannel,massDS.Label);
for i=1:2:coreNum
    dataNoAgg=MibiRemoveAggregates(p{i}.countsNoNoise(:,:,currChannelInd),gausRad,t,gausFlag);
end

