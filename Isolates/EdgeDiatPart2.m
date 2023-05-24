function [radii,centers,biggest] = EdgeDiatPart2(nameSort)

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

nameCharFirst = char(nameSort{1}); % read the first time point
Im1=imread(nameCharFirst,2);
I1 = (Im1);
% ---  Enhancement ---
Ienh =  imadjust(I1);
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
end