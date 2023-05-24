function [BwN]= DiaTindexingID_19(nameSort)

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

% Function purpose
% Sort by multipoint and fluorescence

% Long Chain order of channels
% CH 1 = BF; CH 2 = CY5; CH 3 = ds-RED (a.k.a. TRITC); CH 4 = GFP (a.k.a. FITC); 5 = CFP

% Short Chain order of channels
% CH 1 = BF; CH 2 = ds-RED (a.k.a. TRITC); CH 3 = GFP (a.k.a. FITC);

% for k = 1:length(nameSort)
thisFileName = nameSort{end}; % Get the full or base file name somehow.
if startsWith(thisFileName, '._')
else
    nameChar = char(thisFileName);
    I_nar=imread(nameChar,1); % read the image
    J = imadjust(I_nar);
    I_narblur = imgaussfilt(J,1);
    [~,threshold] = edge(I_narblur,'sobel');
    fudgeFactor = 0.8; %0.9
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
    BwN = BWfinal;
end
end



