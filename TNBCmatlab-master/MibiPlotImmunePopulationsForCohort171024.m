% MibiPlotImmunePopulationsForCohort170911
% plot classification of cells to immune groups for all patients

massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg';
pathSegment = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/SegmentDeep';
pathCluster = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResultsDeep/';
imSize=2048;
channelNum = length(massDS);
points=43;
load([pathSegment,'/dataWithTumorAndImmuneGroups171024.mat']);
T=50;
tumorYNcol = 53;
groupCol=57;


images=zeros(imSize,imSize,3,points,'uint8');

% colorcode:
% [ 0 bg - white, 0 tumor  - gray, immune- rainbow ]
cmap = [200 200 200; 127 0 255; 255 0 255; 255 0 127; 255 0 0; 255 127 0; 255 255 0; 127 255 0; 0 255 0; 0 255 127; 0 255 255; 0 127 255; 0 0 255];
cmap01 = cmap/255;

for i=[1:29,31:points+1]
    disp(i);
    load([pathSegment,'/Point',num2str(i),'/segmentationParams.mat']);
    stats = regionprops(newLmod,'PixelIdxList');
    currInds = (dataAll(:,1) == i);
    currCellData = dataAll(currInds,:);
    currIndsLabel = (labelIdentityAll(:,1) == i);
    currLabelIdentity = labelIdentityAll(currIndsLabel,:);
    labelNum = length(currLabelIdentity);
    imageL = zeros(size(newLmod,1)*size(newLmod,2),1);
    counts=0;
    for j=1:labelNum-1
        cellInd = find(ismember(currCellData(:,2),j));
        if ~isempty(cellInd)
            % check if tumor
            if (currCellData(cellInd,tumorYNcol) == 1)
                cellVal =1;
            else
                cellVal = currCellData(cellInd,groupCol) + 1;
            end
            imageL(stats(j).PixelIdxList)=cellVal;
        end
    end
    imageLReshape = reshape(imageL,size(newLmod,1),size(newLmod,2));
    images([1:size(imageLReshape,1)],[1:size(imageLReshape,2)],:,i) = label2rgb(imageLReshape,cmap01,'w');
    % plot single
    f=figure;
    colormap('parula');
    imagesc(label2rgb(imageLReshape,cmap01,'w'));
    plotbrowser on;
    saveas(f,[path,'/Point',num2str(i),'/TIFs3/immuneGroupsDeep171024.png']);
end

% get order to plot
% load([pathSegment,'/mixingDataDivByTotal170907.mat']);
% 
% plot all together
% f=figure;
% [ha, pos] = tight_subplot(6,7,[.01 .01],[.01 .01],[.01 .01]); 
% for j = 1:points
%     axes(ha(j));
%     data = images(:,:,:,j);
%     imshow(data);
%     title(['Sample ',num2str(j)]);
% end
% saveas(f,[path,'/ByChannel3/groupsByImmuneMixing170911.fig']);