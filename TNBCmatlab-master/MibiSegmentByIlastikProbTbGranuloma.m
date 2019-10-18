    %% Pipeline for nuclear segmentation using pixel probabilities from Ilastik


    %% Get perimiters of nuclei
    % load data and get nuclear markers
    nucIm = imread('processed_dsDNA.tif');
    imSize = (size(nucIm,1));
    % apply gaussian filter
    nucIm2_1 = wiener2(nucIm,[3 3]);
    nucIm3 = imgaussfilt(nucIm2_1,1);
    nucIm4 = mat2gray(nucIm3);
    nucIm5 = imadjust(nucIm4);
    % binarize 
    bw1 = imbinarize(nucIm3, 'adaptive');
    SE = strel('disk',4);
    bw=imdilate(bw1,SE);
    % get perimiter
    bw_perim = bwperim(bw);
    % plot nuclei with perim
    maxv=25;
    rgb_image = MibiGetRGBimageFromMat(nucIm,maxv);
    rgb_image_perim= imoverlay(rgb_image , bw_perim , [1 .3 .3]);
    figure; imagesc(rgb_image_perim);

    % %% Get maxima from Ilastik probabilities
    % read nuclear segmentation from Ilastik
    %dataIl = h5read([path,'/ForIlastik/IlastikBatchP1/p_',num2str(pointNumber),'_nuclei_Probabilities_underSegment.h5'],'/exported_data');
    dataIl = h5read(['processed_dsDNA2_Probabilities.h5'],'/exported_data');
    probNuc=reshape(dataIl(1,:,:),imSize,imSize);
    probNuc=probNuc';

    % %% Get maxima from Ilastik probabilities by 3 colors
    % % read nuclear segmentation from Ilastik
    % %dataIl = h5read('Point2/ForIlastik/dsdna_panck_cd45_Probabilities.h5','/exported_data');
    % dataIl = h5read([path,'/Point',num2str(pointNumber),'/ForIlastik/nuc_immune_tumor_Probabilities.h5'],'/exported_data');
    % Labels = {'CD45','panCk','dsDNA_tumor','dsDNA_imune','Bg'};
    % probAll = zeros(imSize,imSize,5);
    % for i=1:5
    %     tmp=reshape(dataIl(i,:,:),imSize,imSize);
    %     tmp=tmp';
    %     probAll(:,:,i)=tmp;
    %     figure; imagesc(probAll(:,:,i));
    % end
    % probNuc=max(probAll(:,:,3),probAll(:,:,4));
    % figure; imagesc(probNuc);
    % 
    % % reduce probability of nuclei in places that have high probablity for
    % % other labels
    % probNuc2=probNuc.*(1-probAll(:,:,1));
    % figure; imagesc(probNuc2);
    % probNuc3=probNuc2.*(1-probAll(:,:,2));
    % figure; imagesc(probNuc3);
    % 
    % probNuc = probNuc3;

    % find local maxima in probability map
    % pks = imregionalmax(probNuc);
    % overlay2 = imoverlay(mat2gray(nucIm3) , bw_perim | pks, [1 .3 .3]);
    maxs = imextendedmax(probNuc,0.05);
    rgb_image_perim_extMax = imoverlay(rgb_image_perim , maxs, [1 0 0]);
    %overlay2 = imoverlay(mat2gray(nucIm3) , bw_perim | maxs, [1 .3 .3]);    
    figure;
    imagesc(rgb_image_perim_extMax);
    %L1 = bwlabel(maxs,4);
    % figure;
    % imshow(label2rgb(L1, @jet, [.5 .5 .5],'shuffle'));
    % title('Before splitting');

    % BW2 = bwmorph(maxs,'open');
    % overlay3 = imoverlay(mat2gray(nucIm3) ,bw_perim | BW2, [1 .3 .3]);
    % figure;
    % imshow(overlay3);

    %% watershed
    [B,L] = bwboundaries(maxs,4,'noholes');
    maxsFix = bw1 & maxs;

    % modify the image so that the background pixels and the extended maxima pixels are forced to be the only local minima in the image.
    Jc = imcomplement(nucIm3);
    I_mod = imimposemin(Jc, ~bw | maxsFix);
    L = watershed(I_mod);
    labeledImage = label2rgb(L);
    % figure;
    % imshow(label2rgb(L,'gray','k','shuffle'));

    %% 1. For each label, decide whether it is of a nucleus/ background
    t=40;
    labelNum = length(unique(L(:)));
    labelIdentity = zeros (labelNum,1);
    labelPixelsPercentInNucleiMask = zeros(labelNum,1);
    labelSize = zeros(labelNum,1);

    for i=1:labelNum
        [r c] = find(L==i);
        labelSize(i) = length(r); 
        labelMask = (L==i);
        labelPixelsNumInNucleiMask = sum(bw(labelMask));
        labelPixelsPercentInNucleiMask(i) = labelPixelsNumInNucleiMask / labelSize(i);
        if (labelPixelsPercentInNucleiMask(i) > 0.7)
            labelIdentity(i) = 1;
        end
    end

    % 2. Merge small regions within the nuclei mask with their neighbours
    keepVec = ones(labelNum,1);
    newL = L;
    for i=1:labelNum
        if (labelIdentity(i) == 1) && (labelSize(i) < t)
            % disp(['Removing label ',num2str(i),'. Size: ',num2str(labelSize(i))]);
            % get neighbour with largest border that is also in nuclear region
            [neighbourLabels , neighbouringRegionSize] = MibiGetNeighbourLabels (newL, i);
            found = 0;
            [neighbouringRegionSizeS , neighbouringRegionSizeSInd] = sort(neighbouringRegionSize,'descend');
            neighbourLabelsS = neighbourLabels(neighbouringRegionSizeSInd);
            maxInd = 1;
            while ~found
                mergeLabelId = neighbourLabelsS(maxInd);
                if (~(mergeLabelId == 0) && (labelIdentity(mergeLabelId) == 1))
                    found = 1;
                else
                    maxInd = maxInd+1;
                end
                if (maxInd >length(neighbourLabelsS)) % reached end of neighbours with no good merging candidate
                    disp (['Warning: no good merging target found for label', num2str(i), '. Keeping it.']);
                    break;
                end
            end
            % update
            if (maxInd <= length(neighbourLabelsS))
                [newL] = MibiMergeLabels (newL, i, mergeLabelId);
                keepVec(i) = 0;
            end
        end
    end

    % Update label numbers to account for deleted labels
    allLabels = [1:labelNum];
    currLabels =  allLabels(keepVec == 1);
    labelIdentityNew = zeros(length(currLabels),1); 
    newLmod = newL;
    for i = 1:length(currLabels)
        newLmod(newLmod == currLabels(i)) = i;
        labelIdentityNew(i) = labelIdentity(currLabels(i));
    end

    cellPerimNewMod= bwperim(newLmod);
    cellPerimNewMod= newLmod;
    cellPerimNewMod(newL>0) = 100;
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    maxv=15;
    rgb_image = MibiGetRGBimageFromMat(double(nucIm),maxv);
    rgb_image_perim= imoverlay(rgb_image , cellPerimNewMod , [1 .3 .3]);
    figure; imagesc(rgb_image_perim);


%     com_cellPerimNewMod = imoverlay(com, cellPerimNewMod, [1 1 1]);
%     figure;
%     imagesc(com_cellPerimNewMod);
%    mkdir([resultsPath,'/Point',num2str(pointNumber)]);
%    save([resultsPath,'/Point',num2str(pointNumber),'/segmentationParams.mat'],'newLmod','cellPerimNewMod','labelIdentityNew');
%    close all;

