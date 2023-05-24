function [output,nmbrD,areaDiat] = SelectCountDiatom(biggest,I_BF,AreD,Elong,AreDmax,fudgeFactor,filt)

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

% This function detect CA diatoms that have an elongated shape
% Elongation number is above 1.25.
% ----Usual parameter----
% Elong = 1.25;
% AreD = 100;
% AreDmax = 800;
% fudgeFactor = 0.9; %0.8
% filt = 2;

nameChar = char(I_BF);
Im = imread(nameChar,1);
J = imadjust(Im);
I_narblur = imgaussfilt(J,filt);
[~,threshold] = edge(I_narblur,'sobel');
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
statsA = regionprops(cc);
threshold = 100000;
removeMask = [statsA.Area]>threshold;
BWfinal(cat(1,cc.PixelIdxList{removeMask})) = false;

BWfinal_ID = BWfinal;
BWfinal_ID(~biggest) = 0;
[Bk,Lk] = bwboundaries(BWfinal_ID,'noholes');
cc2 = bwconncomp(BWfinal_ID);
output = zeros(size(BWfinal_ID));

% --- LABELLED OBJECT ---
% properties of ROI boundaries
stats = regionprops(Lk,'Area','Centroid', 'Perimeter','MajorAxisLength','MinorAxisLength');
info = regionprops(BWfinal_ID,'Boundingbox') ;

% % % % % %------------------------------------------------------
% % % % % %---VISUALIZATION--------------------------------------
% % % % % %------------------------------------------------------
% figure,
% imshowpair(J,BWfinal_ID)
% hold on
% % % % % %------------------------------------------------------
coord_diat = {};
counts = {};
AreaCell = {};
% loop over the boundaries
for kj = 1:length(Bk)
    % obtain (X,Y) boundary coordinates corresponding to label 'k'
    boundary = Bk{kj};
    perimeter = stats(kj).Perimeter;
    % obtain the area calculation corresponding to label 'k'
    area = stats(kj).Area;
    % compute the roundness metric
    metric = (perimeter.^ 2) ./ (4 * pi * area);
    Elongation = sqrt(stats(kj).MajorAxisLength/stats(kj).MinorAxisLength);
    % display the results
    metric_string = sprintf('%2.2f',Elongation);
    % mark objects that are near to be circular
    
    if  (Elongation> Elong && area>AreD && area<AreDmax)%30
        
        counts{kj} = 1;
        [y,x] = ind2sub(cc2.ImageSize, cc2.PixelIdxList{kj});
        output(cc2.PixelIdxList{kj}) = 1;
        coord_diat{end+1} = [x, y];
        boundary = Bk{kj};
        area = stats(kj).Area;
        
        % % % % %------------------------------------------------------
        % % % % %---VISUALIZATION--------------------------------------
        % % % % %------------------------------------------------------
        % plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
        % % % % %------------------------------------------------------
        
        % Boundingbox
        BBk = info(kj).BoundingBox;
        
        %         % % % % % %------------------------------------------------------
        %         % % % % % %---VISUALIZATION--------------------------------------
        %         % % % % % %------------------------------------------------------
        %         rectangle('Position', [BBk(1),BBk(2),BBk(3),BBk(4)],'EdgeColor','r','LineWidth',2) ;
        %         text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','r',...
        %             'FontSize',5,'FontWeight','bold')
        %         % % % % % %------------------------------------------------------
        
        AreaCell{kj} = area;
    end
end
areaDiat = AreaCell;
nmbrD = sum([counts{:}]);
end





