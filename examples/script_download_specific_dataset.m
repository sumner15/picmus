clear all;
close all;
clc;


%-- parameters
signal_selection = 2;       %-- 1: RF | 2: IQ
pht_selection = 1;          %-- 1: numerical | 2: in_vitro_type1 | 3: in_vitro_type2 | 4: in_vitro_type3
transmission_selection = 1; %-- 1: regular | 2: dichotomous
nbPW_number = 1;           %-- An odd value between 1 and 75


%-- generate corresponding dataset filename
[filenames] = tools.generate_filenames(signal_selection,pht_selection,transmission_selection,nbPW_number);


%-- data location
url = 'https://www.creatis.insa-lyon.fr/EvaluationPlatform/picmus/dataset/';
local_path = [picmus_path(),'/data/']; % location of example data in this computer


%-----------------------------------------------------
%-----------------------------------------------------
%-- download scanning region file
if (~exist([local_path,filenames.scan],'file'))
     tools.download(filenames.scan, url, local_path);
end

%-- download dataset
if (~exist([local_path,filenames.dataset],'file'))
    tools.download(filenames.dataset, url, local_path);
end


