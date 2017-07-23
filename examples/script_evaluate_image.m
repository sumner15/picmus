clear all;
close all;
clc;


%-- parameters
signal_selection = 2;       %-- 1: RF | 2: IQ
pht_selection = 1;          %-- 1: numerical | 2: in_vitro_type1 | 3: in_vitro_type2 | 4: in_vitro_type3
transmission_selection = 1; %-- 1: regular | 2: dichotomous
nbPW_number = 1;            %-- An odd value between 1 and 75


%-- generate corresponding dataset filename
[filenames] = tools.generate_filenames(signal_selection,pht_selection,transmission_selection,nbPW_number);


%-- set path
path_img = [picmus_path(),'/results/',filenames.image];


%-- Read reconstructed image
image = us_image();
image.read_file(path_img);


% Setup info structure
info = tools.generate_data_info_structure(filenames.pht_name);


% Create picmus metric data object
metrics = us_picmus_metrics();
metrics.image = image;
metrics.scan = image.scan;
metrics.set_data_information(info);
metrics.flagDisplay = 1; % Set this flag to 1 to display intermediate results
metrics.evaluate(); % This method launch all the metrics at once


%-- Display scores
disp('--------------------------------')
disp('Speckle quality')
disp(metrics.scoreSpeckleQuality)
disp('--------------------------------')
disp('Geometric distortion')
disp(metrics.scoreGeometricalDistortion)
disp('--------------------------------')
disp('Intensity linearity')
disp(metrics.scoreLinearIntensity)
disp('--------------------------------')
disp('Full width at half maximum')
disp(metrics.scoreFWHM)
disp('--------------------------------')
disp('Contrast')
disp(metrics.scoreContrast)
disp('--------------------------------')
disp('Axial resolution')
disp(metrics.scoreResolutionAxial)
disp('--------------------------------')
disp('Lateral resolution')
disp(metrics.scoreResolutionLateral)



