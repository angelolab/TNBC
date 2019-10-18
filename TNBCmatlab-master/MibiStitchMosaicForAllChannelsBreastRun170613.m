function MibiStitchMosaicForAllChannelsBreastRun170613
% function MibiStitchMosaicForAllChannelsBreastRun170613(outputDir,xdTopMed,ydTopMed,xdLeftMed,ydLeftMed,xdRightMed,ydRightMed)
% function gets the stitching parameters for the run and creates a tiff
% directory of all channels stitched. Note that this specific run had a
% special handeling of line four of the run, preventing this script from
% being generic

xNumPoint = 20;
yNumPoint = 26;
pointNum=xNumPoint*yNumPoint;
dataSize= 508;
channelNum=49;
outputDir= 'stitchImages';

xdTopMed=6;
ydTopMed=-65;
xdLeftMed=444;
ydLeftMed=11;
xdRightMed=64;
ydRightMed=-11;

fileList=dir('Point1');

for f=3:length(fileList)
    %% load points in serpentine scheme into 3d mat just for plotting
    % load tifs
    allData = zeros(xNumPoint*dataSize,yNumPoint*dataSize);
    allDataCell = cell(xNumPoint,yNumPoint);
    for i=1:xNumPoint
        xloc=xNumPoint-i+1;
        for j=1:yNumPoint
            yloc=j;
            if (mod(i,2) == 0)
                yloc = yNumPoint-j+1;
            end
            currPoint = (i-1)*yNumPoint+j;
            allDataCell{xloc,yloc} = double(imread(['Point',num2str(currPoint),'/',fileList(f).name]));
            allDataCell{xloc,yloc} = allDataCell{xloc,yloc}([3:dataSize+2],[3:dataSize+2]);
            allData([(xloc-1)*dataSize+1 : (xloc-1)*dataSize+dataSize] , [(yloc-1)*dataSize+1 : (yloc-1)*dataSize+dataSize]) = allDataCell{xloc,yloc};
        end
    end


    %% stitch
    allDataStitch = zeros((xNumPoint+2)*dataSize,(yNumPoint+2)*dataSize);
    startPosGlobal = ([(xNumPoint)*dataSize-dataSize/2,dataSize]);

    for i=1:xNumPoint
        xloc=xNumPoint-i+1;
        direction=1; % serpentine direction right. Coregister left
        for j=1:yNumPoint
            yloc=j;
            if (mod(i,2) == 0)
                yloc = yNumPoint-j+1;
                direction=0; % serpentine direction left. Coregister right
            end
            currPoint = (i-1)*yNumPoint+j;
            % get curr data
            currData = allDataCell{xloc,yloc};
            % get curr position
            if (xloc == xNumPoint) && (yloc==1) % first point. No coregistering
                currPos = startPosGlobal;
            elseif (direction==1) && ~(yloc == 1) % registering along right movement
                lastPos = currPos;
                currPos = ([lastPos(1)+ydLeftMed,lastPos(2)+xdLeftMed]);
                currData([1:dataSize-ydLeftMed],[1:dataSize-xdLeftMed])=0;
                if (i>1) && ~(xloc == 4) % not first row and not row 4
                    currData([dataSize+ydTopMed+1:dataSize],[1:dataSize])=0;
                end
            elseif ((direction==0) && (yloc == yNumPoint)) || ((direction==1) && (yloc == 1)) % registration along top movement
                lastPos = currPos;
                if (xloc == 4)
                    currPos = ([lastPos(1)-dataSize,lastPos(2)]); % for row four do not overlap. Probably something happenned during the run
                else
                    currPos = ([lastPos(1)-dataSize-ydTopMed,lastPos(2)+xdTopMed]);
                    if (yloc == yNumPoint)
                        currData([dataSize+ydTopMed+1:dataSize],[1:dataSize-xdTopMed])=0;
                        currData([dataSize+ydTopMed-ydLeftMed+1:dataSize],[1:dataSize-xdTopMed-xdLeftMed])=0;
                    elseif (yloc==1)
                        currData([xdLeftMed+xdTopMed+1:dataSize],[1:dataSize-xdTopMed])=0;
                        currData([dataSize+ydTopMed+ydLeftMed+1:dataSize],[dataSize-xdTopMed+1:dataSize])=0;
                    end
                end
            elseif (direction==0) && ~(yloc == yNumPoint) % registering along left movement
                lastPos = currPos;
                %currPos = ([lastPos(1)+ydRightMed,lastPos(2)-dataSize+xdRightMed]);
                currPos = ([lastPos(1)+ydRightMed,lastPos(2)-xdLeftMed]);
                currData([dataSize+ydTopMed+1:dataSize],[1:dataSize-xdTopMed])=0;
                currData([ydLeftMed+1:dataSize],[xdLeftMed+1:dataSize])=0;
                if ~(yloc == 1)
                    currData([dataSize+ydTopMed-ydLeftMed+1:dataSize],[1:dataSize-xdTopMed-xdLeftMed])=0;
                end
            end
            % get previous data in those positions
            prevData = allDataStitch([currPos(1):currPos(1)+dataSize-1],[currPos(2):currPos(2)+dataSize-1]);
            % zero out pixels in new data that are positive in prev data
            %currData(prevData>0)=0;
            % update
            allDataStitch([currPos(1):currPos(1)+dataSize-1],[currPos(2):currPos(2)+dataSize-1])=prevData+currData;
        end
    end
    
    %% save
    data2write=uint8(allDataStitch);
    imwrite(data2write,[outputDir,fileList(f).name]);
end

