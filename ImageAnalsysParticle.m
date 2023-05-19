%% ---- Analyse particles with isolates-----

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

clc;
clear;
close all;

% Choose which dataset to analyze

%Folder structure contain multiple XY position and subsequent Z

for brkn = 1:2
    if brkn ==1
        rootDir = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
        rootDirAnaerobic = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
        load 'OxygenCalibration2.mat'
        maxY = median([StructOxygenCalibration.AerobicCal{:}]);
        stdMaxY = std([StructOxygenCalibration.AerobicCal{:}]);
        minY = median([StructOxygenCalibration.AnaerobicCal{:}]);
        stdminY = std([StructOxygenCalibration.AnaerobicCal{:}]);
        OxuM_l = 220.22;
        fluoNANO = [maxY:1:minY];
        oxConc = OxuM_l*((fluoNANO)-minY)./(maxY-minY)
        OxZone = [40, 20, 10, 6.25, 5, 1]; % Oxygen conc. cut off
        ConvDenFluoZones = {}
        for jjh = 1:length(OxZone)
            ConvDenFluoZones{jjh} = min(fluoNANO(oxConc<OxZone(jjh)))
        end
    else
        rootDir = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
        rootDirAnaerobic = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
        
        load 'OxygenCalibration2.mat'
        maxY = median([StructOxygenCalibration.AerobicCal{:}]);
        stdMaxY = std([StructOxygenCalibration.AerobicCal{:}]);
        minY = median([StructOxygenCalibration.AnaerobicCal{:}]);
        stdminY = std([StructOxygenCalibration.AnaerobicCal{:}]);
        OxuM_l = 220.22;
        fluoNANO = [maxY:1:minY];
        oxConc = OxuM_l*((fluoNANO)-minY)./(maxY-minY)
        OxZone = [40, 20, 10, 6.25, 5, 1];
        ConvDenFluoZones = {}
        for jjh = 1:length(OxZone)
            ConvDenFluoZones{jjh} = min(fluoNANO(oxConc<OxZone(jjh)));
        end
        
    end
    
    % ----- ANALISE IMAGES: GROWTH AND LOCAL  OXYGEN ------
    
    StructMultipleROI_ID = struct;
    for i=1:length(rootDir)
        cd(rootDir{i});
        filCorrupted1 = dir('*._T*');
        nameCorrupted1 = {filCorrupted1(:).name};
        for kj = 1:length(nameCorrupted1)
            delete(char(nameCorrupted1{kj}));
        end
        
        %---Define path---
        file=dir([rootDir{i}, '/', '*.tif']);
        name = {file(:).name};
        %nameSort=natsort(name);
        NBrf = zeros(1,length(name));
        nameSort=natsort(name);
        
        rootDirNth = rootDir{i};
        [BwN]= DiaTindexingID_19(nameSort);
        cd(rootDir{i});
        [radii,centers,biggest] = EdgeDiatPart2(nameSort);
        
        % ---- OXNANO AEROBIC ----
        nameCharOX = char(nameSort{1});
        OxImAE = imread(nameCharOX,2); %to be further processed
        OxImAE_Orig = OxImAE; %original no modification
        
        % ---- OXNANO ANAEROBIC ----
        cd(rootDirAnaerobic{i});
        fileAN=dir([rootDirAnaerobic{i}, '/', '*.tif']);
        nameAN = {fileAN(:).name};
        nameCharOXAN = char(nameAN{1});
        OxImAN = imread(nameCharOXAN,2);
        OxImAN_Orig = OxImAN;
        [radiiAN,centersAN,biggestAN] = EdgeDiatPart2(nameAN);
        cd(rootDir{i});
        
        %--------IDENTIFY DIATOM-----
        Elong = 1.25;
        AreD = 100;
        AreDmax = 800;
        fudgeFactor = 0.9;
        filt = 2;
        I_BF = nameSort{1};
        [output,nmbrD] = SelectCountDiatom(biggest,I_BF,AreD,Elong,AreDmax,fudgeFactor,filt);
        outputInit = output;
        
        % --- IDENTIFY ALL THE COLONIES ---
        % "Backtracking" function
        UpperThreshold = 7;
        LowerThreshold = 5;
        [Rect_v,Img,AreaCh_v,DistEdgCh,BWcolonies] = DiatomBoundingBox_ID_19(nameSort,UpperThreshold,LowerThreshold,radii,centers,biggest);
        
        for kk = 1
            growth = {};
            EdgeD = {};
            OxFrCol= {};
            OxIntColony = {};
            OxOutColony = {};
            OxParticleCol = {};
            oxParticle = {};
            
            for j= 1:length(nameSort)
                
                % --- OXNANO IMAGE ----
                nameCharOX = char(nameSort{j});
                OxIm = imread(nameCharOX,2);
                OxIm2 = OxIm;
                
                % --- OXNANO IMAGE ----
                nameCharOX = char(nameSort{1});
                %[medianOx] = OxygenZone2(nameCharOX);
                
                % --- BF IMAGE ----
                nameChar = char(nameSort{j});
                Im=imread(nameChar,1);
                
                % --- GET INDIVIDUAL COLONIES ---
                [BWfinal] = DiatomROI_19(Im);
                
                %--------IDENTIFY DIATOM-----
                % identify diatoms in the oxnano image.
                Elong = 1.25;
                AreD = 100;
                AreDmax = 800;
                fudgeFactor = 0.9; %0.8
                filt = 2;
                I_BF = nameSort{j};
                [output,nmbrD] = SelectCountDiatom(biggest,I_BF,AreD,Elong,AreDmax,fudgeFactor,filt);
                
                % Substract the diatoms from the Oxnano image
                output2 = (~output);
                OxIm2(~output2) = NaN;
                
                % Substract diatom from the first image
                outputInit2 = (~outputInit);
                OxIm2(~outputInit2) = NaN;
                OxIm3 = OxIm2; %get signal of all colonies
                OxIm3(~BWfinal)= NaN;
                
                % Substract the diatoms from the ANAEROBIC Oxnano image
                OxImAN(~output2) = NaN;
                % Substract diatom from the first image
                OxImAN(~outputInit2) = NaN;
                
                % Substract the diatoms from the AEROBIC Oxnano image
                OxImAE(~output2) = NaN;
                % Substract diatom from the first image
                OxImAE(~outputInit2) = NaN;
                
                % Oxygen at particle level, without diatom
                BWfinal2 = (~BWfinal); %invert the binary image of the colony
                OxIm2(~BWfinal2) = NaN; %eliminate colonies from the oxnano image
                OxIm2(~biggest)= NaN;
                OxIm2(OxIm2==0) = NaN;
                
                %---------OX. PARTICLE LEVEL---------
                % Measure oxnano median in 6 zones
                for ijx = 1:6
                    radY = [1:500:radii];
                    imageSize2 = size(OxIm2);
                    ci2 = [centers, radii-radY(ijx)];     % center and radius of circle ([c_row, c_col, r])
                    [xx2,yy2] = ndgrid((1:imageSize2(1))-ci2(1),(1:imageSize2(2))-ci2(2));
                    mask2 = uint8((xx2.^2 + yy2.^2)<ci2(3)^2);
                    %--Mask of zone for calibration anaerobic particle---
                    radYAN = [1:500:radiiAN];
                    imageSize3 = size(OxImAN);
                    ci3 = [centersAN, radiiAN-radYAN(ijx)];     % center and radius of circle ([c_row, c_col, r])
                    [xx3,yy3] = ndgrid((1:imageSize3(1))-ci3(1),(1:imageSize3(2))-ci3(2));
                    mask3 = uint8((xx3.^2 + yy3.^2)<ci3(3)^2);
                    % Signal at particle level
                    croImage2 = OxIm2;
                    croImage2(~mask2)= NaN;
                    % Signal at colony level
                    OxIm3Col = OxIm3;
                    OxIm3Col(~mask2)= NaN;
                    % Signal at particle level AEROBIC CALIBRATION
                    croImage3 = OxImAE;
                    croImage3(~mask2)= NaN;
                    % Signal at particle level ANAEROBIC CALIBRATION
                    croImage4 = OxImAN;
                    croImage4(~mask3)= NaN;
                    % Anaerobic
                    ImVc3 = [croImage4(:)];
                    ImVc3 = ImVc3(~isnan(ImVc3));
                    vectIm3 =(ImVc3(ImVc3>0));
                    minYlocal = nanmedian(vectIm3);
                    % Aerobic
                    ImVc4 = [croImage3(:)];
                    ImVc4 = ImVc4(~isnan(ImVc4));
                    vectIm4 =(ImVc4(ImVc4>0));
                    maxYlocal = nanmedian(vectIm4);

                    %--- Store local calibration ----
                    zoneScal{ijx} = {[minYlocal maxYlocal]};
                    [cnt1,bns1z] = imhist([OxIm3Col(OxIm3Col>0)]); %Signal at colony level
                    bns1= OxuM_l*((bns1z)-double(minYlocal))./(double(maxYlocal)-double(minYlocal));
                    [cnt2,bns2z] = imhist([croImage2(croImage2>0)]); %Signal at particle level
                    bns2= OxuM_l*((bns2z)-double(minYlocal))./(double(maxYlocal)-double(minYlocal));
                    prtOX = nanmedian(double(OxIm2(OxIm2>0)));
                    A = (OxIm2(OxIm2>0));
                    B = A(~isnan(A));
                    prtOXiqr = iqr(double(B));
                    % Oxnano intensity distribution at colony and at
                    % particle level at different zones oxnano median and
                    % interquantile at particle level
                    oxParticle{ijx}= [{[cnt1,bns1]} {[cnt2,bns2]} prtOX prtOXiqr];
                end
                for ii = 1:1:length(Rect_v{1,kk})
                    % get the distance of ROI from the edge of particle
                    xCentroid = Rect_v{1,kk}{1,ii}(1) + Rect_v{1,kk}{1,ii}(3)/2;
                    yCentroid = Rect_v{1,kk}{1,ii}(2) + Rect_v{1,kk}{1,ii}(4)/2;
                    centerRect = [xCentroid yCentroid];
                    EdgeDist = (radii-(pdist2(centers, centerRect)));
                    %--GET THE INDEX AREA OF LOCAL CALIBRATION---
                    radYAN = [1:500:radiiAN];
                    [valan,idxan]=min(abs(radYAN-EdgeDist));
                    %--- Crop the ROI
                    I_bub = imcrop(BWfinal,[Rect_v{1,kk}{1,ii}]);
                    BWnobord = I_bub;
                    %--- Crop the ROI Oxnano
                    I_OnIm = imcrop(OxIm,[Rect_v{1,kk}{1,ii}]);
                    %--------GET LOCAL OXYGEN WITH REAL OX VALUE-----
                    [OxIntensityNC,OxFrac] = LocalOxigen3(ConvDenFluoZones,BWnobord,I_OnIm);
                    minYlocalC =zoneScal{1, idxan}{1, 1}(1);
                    maxYlocalC = zoneScal{1, idxan}{1, 1}(2);
                    %--------STORE COLONY REAL OX VALUE-----
                    OxIntensity= OxuM_l*((OxIntensityNC)-double(minYlocalC))./(double(maxYlocalC)-double(minYlocalC));
                    % ---Measure Area of ROI
                    sb = regionprops(BWnobord,'Area');
                    % ---Conditional to store zero value ---
                    if isempty(sb)
                        AreaB = 0;
                        OxIntensity = 0;
                    else
                        for k = 1:length(sb)
                            AreaB(k) = sb(k).Area;
                        end
                    end
                    
                    if (AreaB==0)
                        maxA = 0;
                    else
                        maxA = max(AreaB);
                    end
                    % % ---- STORE ----
                    maxArea = (maxA);
                    OxFracColony{j,ii} = OxFrac;
                    OxIntColony{j,ii} = OxIntensity;
                    growth{j,ii} = maxArea;
                    EdgeD{j,ii} = EdgeDist;
                end
                OxParticleCol{j} = oxParticle;
            end
            maxArea = (maxA);
            OxFrCol{kk} = OxFracColony;
            OxIntCol{kk} = OxIntColony;
            GrowthID{kk}= growth;
            EdgeDID{kk} = EdgeD;
        end
        StructMultipleROI_ID.ParticleOX = OxParticleCol; %Oxnano particle
        StructMultipleROI_ID.GrowthTime = GrowthID; %Growth
        StructMultipleROI_ID.EdgeDistances = EdgeDID; %Distance to edge
        StructMultipleROI_ID.AreaChannel = AreaCh_v; %Area channel
        StructMultipleROI_ID.DistChannel = DistEdgCh; %Distance to edge channel
        StructMultipleROI_ID.DynamicOx = OxIntCol; %Ox at colony level
        StructMultipleROI_ID.DynamicFractOx = OxFrCol; %Colony oxnano zones
        save('MultipleROI_ID7.mat', 'StructMultipleROI_ID');
    end
end


