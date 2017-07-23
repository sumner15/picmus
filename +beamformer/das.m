function das(scan_path,dataset_path,image_path)

%-- Implementation of conventional Delay-and-Sum (DAS) with apodization in reception
%-- The code works both for IQ and RF format

%-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

%-- $Last Update: 2017/06/28 $

%-- load scan and dataset
scan = linear_scan();
scan.read_file(scan_path);
dataset = us_dataset();
dataset.read_file(dataset_path);

%-- check if data is RF
if(isempty(dataset.modulation_frequency)||dataset.modulation_frequency==0)
    dataset.data=reshape(hilbert(dataset.data(:,:)),size(dataset.data));
end

%-- receive apodization:
%-- dynamically expanding receive aperture with Tukey 25% apodization
rx_f_number = 1.75;
rx_aperture = scan.z/rx_f_number;
rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey25');
receive_apodization = bsxfun(@rdivide,receive_apodization,sum(receive_apodization,2));

%-- angular apodization -> no apodization
angular_apodization = ones(scan.pixels,dataset.firings);
angular_apodization  = bsxfun(@rdivide,angular_apodization,sum(angular_apodization,2));

%-- beamforming loop
beamformed_data = zeros(scan.pixels,1);
time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
w0 = 2*pi*dataset.modulation_frequency;
wb = waitbar(0,'DAS beamforming');

for pw=1:length(dataset.angles)
    %-- transmit delay
    transmit_delay = scan.z*cos(dataset.angles(pw))+scan.x*sin(dataset.angles(pw));
    for nrx=1:dataset.channels
        %-- progress bar
        step=(nrx + (pw-1)*dataset.channels)/length(dataset.angles)/dataset.channels;
        waitbar(step,wb,sprintf('DAS-RF beamforming %0.0f%%',step*100));
        %-- receive delay
        receive_delay = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
        %-- total delay
        delay = (transmit_delay+receive_delay)/dataset.c0;
        %-- phase shift
        phase_shift = exp(1i.*w0*(delay-2*scan.z/dataset.c0));
        %-- beamformed data
        beamformed_data = beamformed_data + phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
    end
end
close(wb);

%-- reshape
envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)]));

%-- declare an us_image object to store the beamformed data
image = us_image('DAS beamforming');
image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
image.algorithm = 'Delay-and-Sum (IQ version)';
image.scan = scan;
image.number_plane_waves = length(dataset.angles);
image.data = envelope_beamformed_data;
image.transmit_f_number = 0;
image.receive_f_number = rx_f_number;
image.transmit_apodization_window = 'none';
image.receive_apodization_window = 'Tukey 25%';
image.write_file(image_path);

end
