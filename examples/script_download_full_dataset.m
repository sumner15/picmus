clear all;
close all;
clc;

%-- data location
url = 'https://www.creatis.insa-lyon.fr/EvaluationPlatform/picmus/dataset/';
local_path = [picmus_path(),'/data/']; % location of example data in this computer


%-----------------------------------------------------
%-----------------------------------------------------
%-- download scanning region file
filename_scan = 'scanning_region_picmus.hdf5';
if (~exist([local_path,filename_scan],'file'))
     tools.download(filename_scan, url, local_path);
end


%-----------------------------------------------------
%-----------------------------------------------------
%-- download database

%-- phantom type
pht_name = cell(1,4);
pht_name{1} = 'numerical';
pht_name{2} = 'in_vitro_type1';
pht_name{3} = 'in_vitro_type2';
pht_name{4} = 'in_vitro_type3';

%-- signal type
signal_format = cell(1,2);
signal_format{1} = 'rf';
signal_format{2} = 'iq';

%-- transmission scheme
transmission_scheme = cell(1,2);
transmission_scheme{1} = '1';
transmission_scheme{2} = '2';

%-- number of steered plane-waves
nbPW = 1:2:75;

%-- Download the full database
step = 0;
total = length(pht_name)*length(signal_format)*length(transmission_scheme)*length(nbPW);
wb = waitbar(0,'Downloading full database');

for i=1:length(pht_name)
    for j=1:length(signal_format)
        for k=1:length(transmission_scheme)
            for l=1:length(nbPW)
                %-- progress bar                
                step = step+1;
                waitbar((step/total),wb,sprintf('Downloading full database %0.0f%%',(step/total)*100));
                %-- set current filename
                filename_dataset = ['dataset_',signal_format{j},'_',pht_name{i},...
                   '_transmission_',transmission_scheme{k},...
                   '_nbPW_',num2str(nbPW(l)),'.hdf5'];
                %-- download file if it has not been done yet
                if (~exist([local_path,filename_dataset],'file'))
                    tools.download(filename_dataset, url, local_path);
                end
           end
       end
   end
end
close(wb);

