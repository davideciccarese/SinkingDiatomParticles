function [Rect_v,Img,AreaCh_v,DistEdgCh,BWcolonies] = DiatomBoundingBox_ID_19(nameSort,UpperThreshold,LowerThreshold,radii,centers,biggest,onOff)


% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

% BACK TRACKING: track individual colonies
nameChar = char(nameSort{end});
Im=imread(nameChar,1); % read the image
I_nar=imread(nameChar,1); % read the image
J = imadjust(I_nar);
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
BWfinal_coord = BWfinal.*biggest;
channels = [1];

Rect_v = {};
AreaCh_v = {};
DistEdgCh = {};

for jz = 1
    
    
    Rect = {};
    AreaCh = {};
    EdgDCh = {};
    
    %FIND ID
    BWfinal_ID = (BWfinal_coord == channels(jz));
    
    % --- LABELLED OBJECT ---
    [Bk,Lk] = bwboundaries(BWfinal_ID,'noholes');
    % properties of ROI boundaries
    stats = regionprops(Lk,'Area','Centroid', 'Perimeter');
    % Boundingbox
    info = regionprops(BWfinal_ID,'Boundingbox') ;
    
    if onOff==1
        % % % % %------------------------------------------------------
        % % % % %---VISUALIZATION--------------------------------------
        % % % % %------------------------------------------------------
        figure,
        %imshowpair(BWfinal_ID,Im)
        imshow(BWfinal_ID) %BWfinal_ID
        hold on
        % % % % %------------------------------------------------------
    else
    end
    
    % loop over the boundaries
    for kj = 321%1:length(Bk)
        % obtain (X,Y) boundary coordinates corresponding to label 'k'
        boundary = Bk{kj};
        perimeter = stats(kj).Perimeter;
        % obtain the area calculation corresponding to label 'k'
        area = stats(kj).Area;
        % compute the roundness metric
        metric = (perimeter.^ 2) ./ (4 * pi * area);
        % display the results
        metric_string = sprintf('%2.2f',metric);
        % mark objects that are near to be circular
        if  ((metric < LowerThreshold) && (metric < UpperThreshold) && area>50)%30 %<------Area is iper critica!
            boundary = Bk{kj};
            area = stats(kj).Area;
            if onOff==1
                % % % % %------------------------------------------------------
                % % % % %---VISUALIZATION--------------------------------------
                % % % % %------------------------------------------------------
                plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
                % % % % %------------------------------------------------------
                hold on
            else
            end
            %
            % Boundingbox
            BBk = info(kj).BoundingBox;
            %
            %             % % % % % %------------------------------------------------------
            %             % % % % % %---VISUALIZATION--------------------------------------
            %             % % % % % %------------------------------------------------------
            %             rectangle('Position', [BBk(1),BBk(2),BBk(3),BBk(4)],'EdgeColor','r','LineWidth',2) ;
            %             %text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','r',...
            %             % 'FontSize',10,'FontWeight','bold')
            %             % % % % % %------------------------------------------------------
            
            if onOff==1
                % % % % % %------------------------------------------------------
                % % % % % %---VISUALIZATION--------------------------------------
                % % % % % %------------------------------------------------------
                hold on
                rectangle('Position', [BBk(1)-20,BBk(2)-20,BBk(3)+20,BBk(4)+20],'EdgeColor','r','LineWidth',2) ;
                text(boundary(1,2)-35,boundary(1,1)+13,sprintf('%2.2f',kj),'Color','r','FontSize',10,'FontWeight','bold')
                % % % % % %------------------------------------------------------
            else
            end
            % get the distance of ROI from the edge of particle
            xCentroid = BBk(1) + BBk(3)/2;
            yCentroid = BBk(2) + BBk(4)/2;
            centerRect = [xCentroid yCentroid];
            EdgeDistCh = (radii-(pdist2(centers, centerRect)));
            Rect{1,kj} = [BBk(1)-20,BBk(2)-20,BBk(3)+20,BBk(4)+20];
            AreaCh{1,kj} = area;
            EdgDCh{1,kj} = EdgeDistCh;
        end
    end
    BWcolonies = BWfinal;
    Img = Im;
    Rect = Rect(~cellfun('isempty',Rect));
    Rect_v{1,jz} = Rect;
    AreaCh_v{1,jz} = AreaCh(~cellfun('isempty',AreaCh));
    DistEdgCh{1,jz} = EdgDCh(~cellfun('isempty',EdgDCh));
end
end
