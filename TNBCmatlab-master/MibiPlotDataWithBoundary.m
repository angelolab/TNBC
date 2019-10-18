function MibiPlotDataWithBoundary(data,boundary,maxVal)
% function MibiPlotDataWithBoundary(data,boundary,maxVal)
% function plots a data matrix, capped at maxVal, overlayed with a boundary

if ~(exist('maxVal')) || isempty(maxVal)
    maxVal = 5;
end
data_rgb = MibiGetRGBimageFromMat(data,maxVal);
data_rgb_perim = imoverlay(data_rgb, boundary, [1 0 0]);
figure; imagesc(data_rgb_perim)