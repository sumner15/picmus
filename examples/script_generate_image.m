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

%-- data location
url = 'https://www.creatis.insa-lyon.fr/EvaluationPlatform/picmus/dataset/';
local_path = [picmus_path(),'/data/']; % location of example data in this computer

%-- download scanning region file
if (~exist([local_path,filenames.scan],'file'))
     tools.download(filenames.scan, url, local_path);
end

%-- download dataset
if (~exist([local_path,filenames.dataset],'file'))
    tools.download(filenames.dataset, url, local_path);
end

%-- set paths
path_dataset = [picmus_path(),'/data/',filenames.dataset];
path_scan = [picmus_path(),'/data/',filenames.scan];
path_image = [picmus_path(),'/results/',filenames.image];


%-- reconstruct image
fprintf(1,'Launching beamforming algorithm...........'); tic;
beamformer.das(path_scan,path_dataset,path_image);
fprintf(1,'done in %0.2fs\n',toc);


%-- visualize beamformed image
image = us_image();
image.read_file(path_image);
image.show(dynamic_range);


disp('Done')

