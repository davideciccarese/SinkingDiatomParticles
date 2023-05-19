%% ---- Aerobic and Anaerobic calibration-----

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

clc;
clear;
close all;

% Choose which dataset to analyze

%brkn = 1; %Intact diatom
%brkn = 2 %Broken diatom

StructOxygenCalibration = struct;

% Folders structure containing multiple XY positions and subsequent Z

if brkn ==1
    % Folders for intact diatom
    rootDirAerobic = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
    rootDirAnaerobic = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
else
    % Folders for broken diatom
    rootDirAerobic = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
    rootDirAnaerobic = {'/tiff_XY_05/tif_XY_05_Z1/'}; %example folder structure
end

AnaerobicOX = {};
AerobicOX = {};
AnaerobicOXprof = {};
AerobicOXprof ={};
EndPointOx = {};

%---Aaerobic----
for i= 1:length(rootDirAerobic)
    cd(rootDirAerobic{i});
    file=dir([rootDirAerobic{i}, '/', '*.tif']);
    name = {file(:).name};
    nameSort=natsort(name);
    
    % --- OXNANO IMAGE ----
    nameCharOX = char(nameSort{1});
    [medianOx] = OxygenZone2(nameCharOX);
    AerobicOX{i,1}= medianOx
    
end

%---Anaerobic----
for i= 1:length(rootDirAnaerobic)
    cd(rootDirAnaerobic{i});
    file=dir([rootDirAnaerobic{i}, '/', '*.tif']);
    name = {file(:).name};
    nameSort=natsort(name);
    
    for j= 1:length(nameSort)
        % --- OXNANO IMAGE ----
        nameCharOX = char(nameSort{j});
        [medianOx] = OxygenZone2(nameCharOX);
        AnaerobicOX{i,j}= medianOx
    end
end

MedAeCalib = AerobicOX;
MedAnCalib = AnaerobicOX;
StructOxygenCalibration.AnaerobicCal = MedAnCalib;
StructOxygenCalibration.AerobicCal = MedAeCalib;
save('OxygenCalibration2.mat', 'StructOxygenCalibration');


