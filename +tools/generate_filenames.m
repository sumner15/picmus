function [filenames] = ...
    generate_filenames(signal_selection,pht_selection,transmission_selection,nbPW_number)
    

    %-- phantom type selection
    switch pht_selection
        case 1
            pht_name = 'numerical';
        case 2
            pht_name = 'in_vitro_type1';
        case 3
            pht_name = 'in_vitro_type2';
        case 4
            pht_name = 'in_vitro_type3';
    end

    %-- signal format selection
    switch signal_selection
        case 1
            signal_format = 'rf';
        case 2
            signal_format = 'iq';
    end

    %-- transmission scheme selection
    switch transmission_selection
        case 1
            transmission_scheme = '1';
        case 2
            transmission_scheme = '2';
    end


    %-- number of steered plane-waves used for reconstruction selection
    if (rem(nbPW_number,2)==0) nbPW_number=nbPW_number+1; end
    if (nbPW_number>75) nbPW_number=75; end
    if (nbPW_number<1) nbPW_number=1; end
    nbPW = num2str(nbPW_number);


    %-- compose corresponding filename_dataset
    filename_dataset = ['dataset_',signal_format,'_',pht_name,...
       '_transmission_',transmission_scheme,...
       '_nbPW_',nbPW,'.hdf5'];
   
    %-- compose corresponding filename_image
    filename_image = ['image_',signal_format,'_',pht_name,...
        '_transmission_',transmission_scheme,...
        '_nbPW_',nbPW,'.hdf5'];

    %-- create corresponding filename_scan
    filename_scan = 'scanning_region_picmus.hdf5';

    %-- create final structure
    filenames = struct('dataset',filename_dataset, ...
        'image',filename_image, ...
        'scan',filename_scan, ...
        'pht_name',pht_name);
    
end

