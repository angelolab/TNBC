% MibiCleanNuclearChannel
% Clean double nuclei for the cohort

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
channelNum = length(massDS);
points=41;
p=cell(1,points);

for pointNumber=[1,41]
    disp(pointNumber);
    p{pointNumber}=load([path,'/Point',num2str(pointNumber),'/dataDeNoiseCohortNoAgg.mat']);
end

data = p{41}.countsNoNoiseNoAgg(:,:,7);
C = normxcorr2(data,data);

imsize =2048;
offset=50;
figure, imagesc(C([imsize-offset:imsize+offset],[imsize-offset:imsize+offset]));