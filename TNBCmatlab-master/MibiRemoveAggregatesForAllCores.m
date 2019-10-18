corePath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/';
destinationPath = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/';
coreNum = 48;
gausRad = 1;
plotFlag = 0;

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);


for i=44:coreNum
    i
    load([corePath,'/Point',num2str(i),'/dataDeNoiseCohort.mat']);
    countsNoNoiseNoAgg = countsNoNoise;
    for j=7:length(massDS)-2
        gausFlag = massDS.GausFlag(j);
        t = massDS.AggFilter(j);
        countsNoNoiseNoAgg(:,:,j) = MibiRemoveAggregates(countsNoNoise(:,:,j),gausRad,t,gausFlag,plotFlag);
    end
    mkdir([destinationPath,'/Point',num2str(i)]);
    save([destinationPath,'/Point',num2str(i),'/dataDeNoiseCohortNoAgg.mat'],'countsNoNoiseNoAgg');
end

