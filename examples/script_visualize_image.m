clear all;
close all;
clc;


%-- parameters
signal_selection = 2;       %-- 1: RF | 2: IQ
pht_selection = 1;          %-- 1: numerical | 2: in_vitro_type1 | 3: in_vitro_type2 | 4: in_vitro_type3
transmission_selection = 1; %-- 1: regular | 2: dichotomous
nbPW_number = 1;            %-- An odd value between 1 and 75
dynamic_range = 60;


%-- generate corresponding dataset filename
[filenames] = tools.generate_filenames(signal_selection,pht_selection,transmission_selection,nbPW_number);


%-- set path
path_img = [picmus_path(),'/results/',filenames.image];


%-- Read reconstructed image
image = us_image();
image.read_file(path_img);
image.show(dynamic_range);
