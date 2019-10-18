% MibiQuantifyImmunecellsPerPatientInTNBCCohort
% For all patients, get the amount of immune cells from different types and
% their distribution

points=41;
massDS = MibiReadMassData(['/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/TA-459_multipleCores2/TNBC_panel_170505_with_bg.csv']);
path = '/Users/lkeren/Box Sync/Leeat_Share/Data/MIBIData/CleanData/NoAgg/ClusteringResults';
load([path,'/ClusterWithGroups170821.mat'],'phenoS');
phenoCol = 92;
[~, patientTF] = ismember(phenoS.textdata,{'patientNum'});
patientIndCol = find(patientTF);
[~, groupTF] = ismember(phenoS.textdata,{'group1'});
groupCol = find(groupTF);

% for each patient get the total number of immune cells and their type
categories={'Tregs','T helper','T cyto','NK','B','Neu','Mac','DC','DC/Mono','Mono/Neu','gdT'};
for i=1:points
    currInds = (phenoS.data(:,patientIndCol)==i);
    currPatientData = phenoS.data(currInds,:);
    for j=1:length(categories)
        currCategoryInds = (currPatientData(:,groupCol)==j);
        countsPerPatient(j,i) = sum(currCategoryInds);
    end
end

totalPerPatient = sum(countsPerPatient,1);

% plot total
figure;
bar(totalPerPatient);

% plot normalized percentages per population
cmap=[127,0,255; ...
    255,0,255; ...
    255,0,127; ...
    255,0,0; ...
    255,127,0; ...
    255,255,0; ...
    127,255,0; ...
    0,255,0; ...
    0,255,127; ...
    0,255,255; ...
    0,127,255; ];
%    0,0,255];
cmap01= cmap./300;

% plot total. Color by population
figure;
bar(countsPerPatient','stacked');
colormap(cmap01);

% sort by total immune cells
[sVals, sInd] = sort(totalPerPatient);
countsSortedByTotal = countsPerPatient(:,sInd);

% plot total. Color by population
totalSort = totalPerPatient(sInd);
countsSortedByTotalNorm = countsSortedByTotal ./ repmat(totalSort,length(categories),1);

figure;
subplot(1,5,[1:2]);
barh(countsSortedByTotalNorm','stacked');
colormap(cmap01);
set(gca,'yTick',[1:points],'yTickLabel',sInd);
ylabel('Patients sorted by total immune infiltrate','fontweight','bold');
axis tight;
title('Composition of immune cells');

subplot(1,5,3);
barh(totalSort,'FaceColor',[0.5,0.5,0.5]);
axis tight;
set(gca,'yTick',[]);
title('Total immune cells');

subplot(1,5,4);
barh(countsSortedByTotalNorm(2,:),'FaceColor',cmap01(2,:));
axis tight;
set(gca,'yTick',[]);
title('% T helper cells');

subplot(1,5,5);
barh(countsSortedByTotalNorm(7,:),'FaceColor',cmap01(7,:));
axis tight;
set(gca,'yTick',[]);
title('% Macrophages');

% get correlation for % of cd4 t-cells and macrophages with total immune
[r1 p1] = corr(totalSort',countsSortedByTotalNorm(2,:)');
[r2 p2] = corr(totalSort',countsSortedByTotalNorm(7,:)');



