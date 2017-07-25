clear all;
close all;
clc;


%-- Constant values corresponding to the probe simulated settings
c0 = 1540;              %-- Speed of sound [m/s]
fs = 100e6;             %-- Sampling frequency [Hz]
f0 = 5.208e6;           %-- Transducer center frequency [Hz]
f_number = 1.75;        %-- f-number used to set the maximum angles
lambda = c0/f0;         %-- Wavelength [m]
pitch = 0.300e-3;       %-- pitch [m]
N_elements = 128;       %-- Number of elements
axialRes = 3 * lambda / 2;
lateralRes = 1.206 * lambda * f_number;


%-- Create a phantom object
pht = us_phantom();
pht.author = 'Olivier Bernard <olivier.bernard@creatis.insa-lyon.fr>';
pht.affiliation = 'CREATIS laboratory - France';


%-- Generate numerical phantom
bckDensity = 20;
padding = 2*1e-3;
probe_width = (N_elements-1) * pitch;
minX = -probe_width/2 - padding;
maxX = probe_width/2 + padding;
minZ = 5*1e-3;
maxZ = 50*1e-3 + padding;
pht.xLimits = [minX maxX];
pht.zLimits = [minZ maxZ];
pht.axialResolution = axialRes;
pht.lateralResolution = lateralRes;
pht.name = 'Picmus calibrating phantom';
pht.setPicmusPhtMode();
pht.bckDensity = bckDensity;
pht_name = 'PicmusPhantom';  
pht.generate();


%-- Save phantom
path_pht = [picmus_path(),'/data/numerical_pht.hdf5'];
pht.write_file(path_pht);


%-- Display phantom
figure; scatter(pht.sca(:,1)*1e3,pht.sca(:,3)*1e3,1,pht.amp(:),'linewidth',3);  
axis equal; grid on; axis ij;
axis([pht.xLimits(1) pht.xLimits(2) pht.zLimits(1) pht.zLimits(2)]*1e3); colorbar;
title('Numercical phantom');   


%-- Perform simulation using FieldII (The folder containing fieldII
%-- input files should be added to MATLAB's path.
sca = pht.sca;
amp = pht.amp;


%-- field II initialisation
field_init(0);
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements


%-- transducer definition P4-2v Verasonics
dt = 1/fs;            %-- Sampling step [s] 
bw = 0.67;            %-- probe bandwidth [1]
height = 5e-3;        %-- Height of element [m]
pitch = 0.300e-3;     %-- pitch [m]
kerf = 0.030e-3;      %-- gap between elements [m]
width = pitch-kerf;   %-- Width of element [m]
elevation_focus = 20e-3;      %-- position of the elevation focus
N_elements = 128;     %-- Number of elements
pulse_duration = 2.5; %-- pulse duration [cycles]
rf_sampling_frequency = 4*f0;          %-- sampling frequency [Hz]
N_plane_waves = 75;                    %-- number of plane waves
alpha_max = 1/2/f_number;
alpha_max = round(alpha_max*180/pi);   %-- to get a positive max angle value in degree
alpha_max = alpha_max*pi/180;          %-- maximum angle [rad]
alpha = linspace(-alpha_max,alpha_max,N_plane_waves).'; %-- vector of angles [rad]


%-- pulse definition
t0 = (-1/bw/f0): dt : (1/bw/f0);
impulse_response = gauspuls(t0, f0, bw);
te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2;   


%-- aperture objects
noSubAz = round(width/(lambda/8));      %-- number of subelements in the azimuth direction
noSubEl = round(height/(lambda/8));     %-- number of subelements in the elevation direction
if exist('elevation_focus','var')
    Th = xdc_focused_array (N_elements, width, height, kerf, elevation_focus, noSubAz, noSubEl, [0 0 Inf]); 
    Rh = xdc_focused_array (N_elements, width, height, kerf, elevation_focus, noSubAz, noSubEl, [0 0 Inf]); 
else
    Th = xdc_linear_array (N_elements, width, height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
    Rh = xdc_linear_array (N_elements, width, height, kerf, noSubAz, noSubEl, [0 0 Inf]);
end


%-- setting excitation, impulse response and baffle
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);


%-- get geometrical center of elements
data = xdc_get(Th,'rect');
geo=data(24:26,:);
x0=zeros(1,N_elements);
for n=1:N_elements
    n_ini=(n-1)*noSubAz*noSubEl+1;
    n_fin=n_ini+noSubAz*noSubEl-1;
    x0(n)=mean(geo(1,n_ini:n_fin));
end

%-- compute CPWC dataset
x_point = sca(:,1);
z_point = sca(:,3);
dt_rf = 1/rf_sampling_frequency;
max_time = 2*sqrt(max([max(x_point)-min(x0) max(x0)-min(x_point)])^2+max(z_point)^2)/c0;      %-- maximum time, samples after this time will be dumped


%-- output data
rf_time_vector = 0:dt_rf:max_time;
CPWC_RF = zeros(length(rf_time_vector),N_elements,N_plane_waves);
band_pass_filter = f0 * (1+bw*[-0.7 -0.5 0.5 0.7]);


%-- compute CPWC signals
for pw=1:N_plane_waves
    % Transmit aperture
    xdc_apodization(Th,0,ones(1,N_elements));
    xdc_times_focus(Th,0,x0(:).'*sin(alpha(pw))/c0);
    % Receive aperture
    xdc_apodization(Rh, 0, ones(1,N_elements));
    xdc_focus_times(Rh, 0, zeros(1,N_elements));    
    % Do calculation
    [v,t] = calc_scat_multi(Th, Rh, sca, amp);
    % Band pass filtering
    v_filtered = tools.band_pass(v,fs,band_pass_filter);
    % lag compensation & decimation
    t_in = (0:dt:((size(v,1)-1)*dt))+t-lag*dt;
    v_compensated = interp1(t_in,v_filtered,rf_time_vector,'pchip',0);    
    CPWC_RF(:,:,pw) = v_compensated;    
    disp(['simulation of rawdata from plane wave n°',num2str(pw),' / ',num2str(N_plane_waves)])
end

%-- normalized the data
CPWC_RF = CPWC_RF / max( CPWC_RF(:) );


%-- create RF dataset objects
rf_ds = us_dataset(pht.name);
rf_ds.probe_geometry = [x0.' zeros(N_elements,2)];      %-- probe geometry [x, y, z] (m)
rf_ds.c0 = c0;                                          %-- refrence speed of sound (m/s)
rf_ds.angles = alpha;                                   %-- angle vector [rad]
rf_ds.data = CPWC_RF;                                   %-- data
rf_ds.initial_time = rf_time_vector(1);                 %-- initial time (s)
rf_ds.sampling_frequency = rf_sampling_frequency;       %-- sampling frequency (Hz)
path_dataset_rf = [picmus_path(),'/data/dataset_rf_numerical.hdf5'];
rf_ds.write_file(path_dataset_rf);


%-- convert RF to IQ
dt_rf = 1/(double(rf_ds.sampling_frequency));
demodulation_frequency = f0;
iq_decimation_factor = round(double(rf_ds.sampling_frequency)/f0);
[CPWC_IQ, iq_time_vector] = tools.demodulate(rf_time_vector.',CPWC_RF,demodulation_frequency,iq_decimation_factor);    


%-- create IQ dataset objects
iq_ds = us_dataset(pht.name);
iq_ds.probe_geometry = [x0.' zeros(N_elements,2)];  %-- probe geometry [x, y, z] (m)
iq_ds.c0 = c0;
iq_ds.angles = alpha;
iq_ds.data = CPWC_IQ;
iq_ds.initial_time = rf_time_vector(1);
iq_ds.sampling_frequency = (double(rf_ds.sampling_frequency)/iq_decimation_factor);
iq_ds.modulation_frequency = demodulation_frequency;
path_dataset_iq = [picmus_path(),'/data/dataset_iq_numerical.hdf5'];
iq_ds.write_file(path_dataset_iq);


%-- Test of reconstruction using standard das method

%-- download scanning region file
url = 'https://www.creatis.insa-lyon.fr/EvaluationPlatform/picmus/dataset/';
local_path = [picmus_path(),'/data/']; % location of example data in this computer
if (~exist([local_path,'scanning_region_picmus.hdf5'],'file'))
     tools.download('scanning_region_picmus.hdf5', url, local_path);
end

%-- set paths
path_dataset = path_dataset_iq;
path_scan = [picmus_path(),'/data/scanning_region_picmus.hdf5'];
path_image = [picmus_path(),'/results/image_iq_numerical.hdf5'];


%-- reconstruct image
fprintf(1,'Launching beamforming algorithm...........'); tic;
beamformer.das(path_scan,path_dataset,path_image);
fprintf(1,'done in %0.2fs\n',toc);


%-- visualize beamformed image
image = us_image();
image.read_file(path_image);
dynamic_range = 60;
image.show(dynamic_range);
disp('Done')

