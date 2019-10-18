% plot all
capImage=5;
%vec=[3,16,17,18,22,36];
%vec=[3,10,40];
vec=[7,40,42,33,21,45,38];

%for i=1:length(massDS)
for j=[1]
    for i=vec
        if i==7
            capImage=25;
        else
            capImage=5;
        end
        data = countsAllSFiltCRSum(:,:,i);
        data(data>capImage)=capImage;
        figure;
        imagesc(data);
        plotbrowser on;
        title(massDS.Label(i));
    end
end
