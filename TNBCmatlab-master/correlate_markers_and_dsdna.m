% plot localization of Ki67,FoxP3,HLA-1,CD45 and dsDNA
%pointNumber=1;
%load(['Point',num2str(pointNumber),'/dataDeNoise.mat']);

marker1List = {'Ki67','FoxP3','CD45','HLA_Class_1'};
marker2List = {'dsDNA'};
[tf1 loc1] = ismember(marker1List,massDS.Label);
[tf2 loc2] = ismember(marker2List,massDS.Label);

rvec = zeros(length(marker1List),1);
pvec = zeros(length(marker1List),1);
distArray = cell(length(marker1List),1);
percentOverlap = zeros(length(marker1List),1);

for i=1:length(marker1List)
    marker1Mat = countsNoNoiseGaus(:,:,loc1(i));
    marker2Mat = countsNoNoiseGaus(:,:,loc2);
%     % correlation
%     % limit the analysis to pixels that are positive to marker 1
%     marker1data = marker1Mat(:);
%     marker2data = marker2Mat(:);
%     marker1dataFilt = marker1data(marker1data>0);
%     marker2dataFilt = marker2data(marker1data>0);
%     [rvec(i) pvec(i)] = corr(marker1dataFilt,marker2dataFilt);
%     
%     % scatter figure
%     figure;
%     scatter(marker1dataFilt,marker2dataFilt);
%     xlabel(marker1List{i});
%     ylabel(marker2List{1});
%     title (['r = ', num2str(rvec(i),2), ' p = ',num2str(pvec(i),2)]);
%     
%     % contour plot
%     % first, get some random x,y coordinates between -3 and 3
%     % to allow the peaks() function to look somewhat meaningful
%     x = rand(10,1)*6 - 3;
%     y = rand(10,1)*6 - 3;
%     % calculate the peaks function for this points
%     z = peaks(x,y);
%     % now, decide the grid we want to see, -3 to 3 at 0.1 intervals
%     % will be fine for this crude test
%     xnodes = -3:0.1:3;
%     ynodes = -3:0.1:3;
%     % now, all gridfit, notice, no sorting!  no nothing!
%     % just tell it the rectangular grid you want, and give it
%     % the raw data.  It will use linear algebra and other robust
%     % techniques to fit the function to the grid points
%     [zg,xg,yg] = gridfit(x,y,z,xnodes,ynodes);
%     % finally, plot the data, Viola!
%     contour(xg,yg,zg)
    
    % distance
    %distArray{i}=MibiGetIntNormDist(countsNoNoise(:,:,loc2),countsNoNoise(:,:,loc1(i)),10,1,10);
    
    % percent overlap
    % limit the analysis to pixels that are positive to marker 1
    marker1data = marker1Mat(:);
    marker2data = marker2Mat(:);
    max(marker1data)
    max(marker2data)
    marker1dataFilt = marker1data(marker1data>0.5);
    marker2dataFilt = marker2data(marker1data>0.5);
    marker2dataFiltPos = marker2dataFilt(marker2dataFilt>10);
    percentOverlap(i) = length(marker2dataFiltPos)/length(marker1dataFilt);
end

%% plot overlap
figure;
bar(percentOverlap);

%% plot distances
% figure;
% for i=1:length(marker1List)
%     histogram(distArray{i},'Normalization','probability','BinEdges',[0:0.1:2],'DisplayStyle','stairs');
%     hold on;
% end
% legend(marker1List);
