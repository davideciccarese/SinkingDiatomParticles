function [BWfinal] = DiatomROI_19(Im)

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

J = imadjust(Im);
I_narblur = imgaussfilt(J,1);
[~,threshold] = edge(I_narblur,'sobel');
fudgeFactor = 0.8;
BWs = edge(I_narblur,'sobel',threshold * fudgeFactor);
se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdil = imdilate(BWs,[se90 se0]);
% --- Fill interior gap ---
BWdfill = imfill(BWsdil,'holes');
% --- Get rid of small pixels ---
seD = strel('diamond',3);
BWfinal = imerode(BWdfill,seD);
% --- Remove big Object ---
% bubble outside of the particle are removed
cc = bwconncomp(BWfinal);
stats = regionprops(cc);
threshold = 100000;
removeMask = [stats.Area]>threshold;
BWfinal(cat(1,cc.PixelIdxList{removeMask})) = false;

end