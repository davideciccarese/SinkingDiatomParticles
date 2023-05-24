function [medianOx] = OxygenZone2(nameCharOX)

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

%This function sample oxnano fluorescence throughout the whole particle

Im1=imread(nameCharOX,2);
I1 = (Im1);
% ---  Enhancement ---
Ienh =  imadjust(I1);
% imshow(Ienh)
% --- Gaussian filter ---
Imgau = imgaussfilt(I1,20);
% --- Otsu threshold ---
Bw = imbinarize(Imgau);
% --- Fill interior gap ---
BWdfill = imfill(Bw,'holes');
% --- ELIMINATE NON CONNECTED ----
biggest = bwareafilt(BWdfill, 1, 'largest');
% --- RADIUS AND CENTER OF THE COLONY ---
stats = regionprops('table',biggest,'Centroid',...
    'MajorAxisLength','MinorAxisLength', 'Area');
% % --- GET THE RADII IN CASE OF MULTIPLE RADII ---
areaMax = max(stats.Area);
C1=stats.Centroid(:,1);
C2=stats.Centroid(:,2);
C=C1(find(stats.Area==areaMax));
CC=C2(find(stats.Area==areaMax));
centers=[C,CC];
majAx= stats.MajorAxisLength(find(stats.Area==areaMax));
minAx= stats.MinorAxisLength(find(stats.Area==areaMax));
diameters = mean([majAx minAx],2);
radii = diameters/2;
% make a binary mask
% convert zeros to NaN
ox1bw = double(biggest);
ox1bw(ox1bw==0) = NaN;
% mask outside the particle with NaN
ox2 = double(I1).*ox1bw;
% %     imagesc(ox2)
% %     imshow(I1)
medianOx = nanmedian([ox2(:)]);
end