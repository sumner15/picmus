classdef us_picmus_metrics < handle
    
    %-- Class defining a procedure to compute the different metrics linked to the PICMUS numerical phantom, i.e.:
    %-- contrast, FWHM, resolution, distortion and linearity of the intensity dynamic,
    %-- the final goal being the assessment of the performance of beamforming techniques in ultrasound imaging
    
    %-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
    %--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

    %-- $Date: 2017/05/23 $      
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class attributes
    
    properties (SetAccess = public) 
        
        %-- administration
        name                         %-- String containing the name of the us_field_simulation
        author                       %-- String containing the name of the author(s) of the beamformed image
        affiliation                  %-- String containing the affiliation of the author(s) of the beamformed image
        creation_date                %-- String containing the date the reconstruction was created

        %-- Common attributes
        pht
        scan
        image
        flagDisplay
        dynamic_range
        scoreSpeckleQuality
        scoreResolutionAxial
        scoreResolutionLateral
        scoreGeometricalDistortion
        scoreLinearIntensity
        scoreFWHM
        scoreContrast
        speckleQualityRoiPsfTimeX
        speckleQualityRoiPsfTimeZ
        speckelQualityRoiCenterX
        speckelQualityRoiCenterZ
        constrastOcclusionCenterX
        constrastOcclusionCenterZ
        constrastOcclusionRadius
        padding
        geometricalDistortionHalfWidthROI1
        geometricalDistortionHalfWidthROI2
        geometricalDistortionX
        geometricalDistortionZ        
        linearIntensityRegionX
        linearIntensityRegionZ
        minX_full
        maxX_full
        fwhmHalfWidthROI
        fwhmX
        fwhmZ
        resolutionAxialPointsX
        resolutionAxialPointsZ
        resolutionLateralPointsX
        resolutionLateralPointsZ
        resolutionHalfWidthROI
        resolutionLeftToRight
        padResolution
        system_lateralResolution
        system_axialResolution
        system_lambda
        system_z_correction
        
    end
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class methods
    
    %-- Constructor
    methods (Access = public)
        
        function h = us_picmus_metrics(input_name)
            if exist('input_name') 
                h.name= input_name; 
            else
                h.name = ' ';
            end            
            h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.author = ' ';
            h.affiliation = ' ';
            h.dynamic_range = 60;
            h.system_z_correction = 0.2e-3;    %-- due to lens curvature
            h.flagDisplay = 0;
            
            %-- Define system/probe settings
            h.system_lambda = 2.9570e-04;   %-- lambda (c/f0)
            h.system_lateralResolution = 0.4;
            h.system_axialResolution = 0.3;
                        
            %-- Define regions to compute speckle quality
            h.speckelQualityRoiCenterX =  [-8  8 -8 10 -8  6 -8  8]*1e-3;
            h.speckelQualityRoiCenterZ =  [12 12 18 20 34 32 40 40]*1e-3;
            h.speckleQualityRoiPsfTimeX = [10 10 10  5 10  5 10 10];
            h.speckleQualityRoiPsfTimeZ = [ 3  3  3  6  3  6  3  3];
            
            %-- Set contrast attributes
            h.padding = 1;
            h.constrastOcclusionCenterX = -8e-3; %-- Position of cyst in x [m]
            h.constrastOcclusionCenterZ = 26e-3; %-- Position of cyst in z [m]
            h.constrastOcclusionRadius = 5e-3;   %-- Radius of cyst [m]
            
            %-- Set geometrical distortion attributes            
            h.geometricalDistortionHalfWidthROI1 = h.system_lambda;            
            h.geometricalDistortionHalfWidthROI2 = 1.8e-3;
            h.geometricalDistortionX = [-15,-7.5,0,7.5,15,-15,-7.5,0,7.5,15]*(1e-3);   %-- x-coordinate of the scatterers [m]
            h.geometricalDistortionZ = [ones(1,5)*(15e-3) ones(1,5)*(37e-3)];          %-- z-coordinate of the scatterers [m]
            
            %-- Set linear intensity attributes
            h.linearIntensityRegionX = [-10 10]*(1e-3);
            h.linearIntensityRegionZ = [45 48]*(1e-3);
            h.minX_full = -19e-3;
            h.maxX_full = +19e-3;

            %-- Set FWHM attributes
            h.fwhmHalfWidthROI = 1.8e-3;
            h.fwhmX = [-15,-7.5,0,7.5,15,-15,-7.5,0,7.5,15]*(1e-3);   %-- x-coordinate of the scatterers [m]
            h.fwhmZ = [ones(1,5)*(15e-3) ones(1,5)*(37e-3)];          %-- z-coordinate of the scatterers [m]

            
            %-- Set resolution attributes
            h.resolutionHalfWidthROI = h.system_lambda/2;
            h.padResolution = h.system_lambda*1.5;
            h.resolutionAxialPointsX = [2.5 3  4  6  9   13   ]*(1e-3);            %-- x-coordinate of the scatterers used to assess the lateral resolution [m]
            h.resolutionAxialPointsZ = [22 23 24 25 25.5 25.75]*(1e-3);            %-- z-coordinate of the scatterers used to assess the lateral resolution [m]                
            h.resolutionLateralPointsX = [2.25 2.5 3 4 6 9 13]*(1e-3);             %-- x-coordinate of the scatterers used to assess the axial resolution [m]
            h.resolutionLateralPointsZ = ones(1,length(h.resolutionAxialPointsX)+1)*(26e-3);    %-- z-coordinate of the scatterers used to assess the axial resolution [m]
            h.resolutionLeftToRight = 1;
            
        end
        
    end
    
    %-- set methods, input format check
    methods  
        
        %-- name
        function set.name(h,input)
            assert(isstr(input), 'Wrong format of the beamformed data name. It should be a string.');
            h.name = input;
        end
        %-- author
        function set.author(h,input)
            assert(isstr(input), 'Wrong format of the author. It should be a string.');
            h.author = input;
        end
        %-- affiliation
        function set.affiliation(h,input)
            assert(isstr(input), 'Wrong format of the affiliation. It should be a string.');
            h.affiliation = input;
        end
        %-- creation_date
        function set.creation_date(h,input_date)
            assert(isstr(input_date), 'Wrong format of the creation date. It should be a string.');
            h.creation_date = input_date;
        end
        %-- pht
        function set.pht(h,input)
            assert(isa(input,'us_phantom'), 'Wrong format of the phantom. It should be a US_PHANTOM class.');
            h.pht = input;
            h.system_lateralResolution = h.pht.lateralResolution;
            h.system_axialResolution = h.pht.axialResolution;
        end
        %-- scan
        function set.scan(h,input)
            assert(isa(input,'linear_scan'), 'Wrong format of the scan. It should be a LINEAR_SCAN class.');
            h.scan = input;
        end
        %-- image
        function set.image(h,input)
            assert(isa(input,'us_image'), 'Wrong format of the image. It should be a US_IMAGE class.');
            h.image = input;
        end
        %-- flag for display
        function set.flagDisplay(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the flag for display. It should be a numeric scalar');
            h.flagDisplay = input;
        end
        %-- flag for display
        function set.dynamic_range(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the dynamic range. It should be a numeric scalar');
            h.dynamic_range = input;
        end
        function set.system_axialResolution(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the axial resolution of the system. It should be a numeric scalar');
            h.system_axialResolution = input;
        end
        function set.system_lateralResolution(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the lateral resolution of the system. It should be a numeric scalar');
            h.system_lateralResolution = input;
        end    
        function set.system_z_correction(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the z correction of the system. It should be a numeric scalar');
            h.system_z_correction = input;
        end        
        function set.constrastOcclusionRadius(h,input_radius)
            assert(size(input_radius,2)==numel(input_radius)&&isnumeric(input_radius), 'Wrong format of the hypoechogenic cyst radius vector. It should be a numeric column vector in (m)');
            h.constrastOcclusionRadius = input_radius;
        end   
        function set.constrastOcclusionCenterX(h,input_centerX)
            assert(size(input_centerX,2)==numel(input_centerX)&&isnumeric(input_centerX), 'Wrong format of the hypoechogenic cyst X center vector. It should be a numeric column vector in (m)');
            h.constrastOcclusionCenterX = input_centerX;
        end 
        function set.constrastOcclusionCenterZ(h,input_centerZ)
            assert(size(input_centerZ,2)==numel(input_centerZ)&&isnumeric(input_centerZ), 'Wrong format of the hypoechogenic cyst Z center vector. It should be a numeric column vector in (m)');
            h.constrastOcclusionCenterZ = input_centerZ;
        end 
        function set.geometricalDistortionX(h,input_distortionX)
            assert(size(input_distortionX,2)==numel(input_distortionX)&&isnumeric(input_distortionX), 'Wrong format of the distortion X center vector. It should be a numeric column vector in (m)');
            h.geometricalDistortionX = input_distortionX;
        end
        function set.geometricalDistortionZ(h,input_distortionZ)
            assert(size(input_distortionZ,2)==numel(input_distortionZ)&&isnumeric(input_distortionZ), 'Wrong format of the distortion Z center vector. It should be a numeric column vector in (m)');
            h.geometricalDistortionZ = input_distortionZ;
        end        
        function set.speckelQualityRoiCenterX(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the speckle quality X center vector. It should be a numeric column vector in (m)');
            h.speckelQualityRoiCenterX = input;
        end
        function set.speckelQualityRoiCenterZ(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the speckle quality Z center vector. It should be a numeric column vector in (m)');
            h.speckelQualityRoiCenterZ = input;
        end        
        function set.speckleQualityRoiPsfTimeX(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the speckle quality X half width vector. It should be a numeric column vector in (m)');
            h.speckleQualityRoiPsfTimeX = input;
        end
        function set.speckleQualityRoiPsfTimeZ(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the speckle quality Z half width vector. It should be a numeric column vector in (m)');
            h.speckleQualityRoiPsfTimeZ = input;
        end        
        function set.fwhmX(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the FWHM X center vector. It should be a numeric column vector in (m)');
            h.fwhmX = input;
        end
        function set.fwhmZ(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the FWHM X center vector. It should be a numeric column vector in (m)');
            h.fwhmZ = input;
        end
        function set.linearIntensityRegionX(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the linear intensity region X center vector. It should be a numeric column vector in (m)');
            h.linearIntensityRegionX = input;
        end
        function set.linearIntensityRegionZ(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the linear intensity region Z center vector. It should be a numeric column vector in (m)');
            h.linearIntensityRegionZ = input;
        end
        function set.resolutionAxialPointsX(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the resolution axial point X vector. It should be a numeric column vector in (m)');
            h.resolutionAxialPointsX = input;
        end
        function set.resolutionAxialPointsZ(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the resolution axial point Z vector. It should be a numeric column vector in (m)');
            h.resolutionAxialPointsZ = input;
        end
        function set.resolutionLateralPointsX(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the resolution lateral point X vector. It should be a numeric column vector in (m)');
            h.resolutionLateralPointsX = input;
        end
        function set.resolutionLateralPointsZ(h,input)
            assert(size(input,2)==numel(input)&&isnumeric(input), 'Wrong format of the resolution lateral point Z vector. It should be a numeric column vector in (m)');
            h.resolutionLateralPointsZ = input;
        end
        function set.resolutionLeftToRight(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the resolution left2right information. It should be a numeric scalar');
            h.resolutionLeftToRight = input;
        end        
        
    end
    
    %-- Public methods
    methods (Access = public)
        

        %-- Generic method called to compute all the metrics
        function evaluate(h)

            %-- Assess speckle quality
            h.evaluateSpeckleQuality();
            
            %-- Assess geometrical distortion
            h.evaluateGeometricalDistortion();
            
            %-- Assess intensity linearity
            h.evaluateLinearIntensity();
            
            %-- Assess contrast
            h.evaluateContrast();
            
            %-- Assess FWHM
            h.evaluateFWHM();
            
            %-- Assess resolution
            h.evaluateResolution();
            
        end


        %------------------------------------------------------
        %-- Assess speckle quality
        function evaluateSpeckleQuality(h)
  
            if (~isempty(h.speckelQualityRoiCenterX))
            
                disp('Perform speckle quality test')
                
                %-- Define parameters / variables
                nb_frames = length(h.image.number_plane_waves);
                frame_list = 1:nb_frames;
                h.scoreSpeckleQuality = zeros(nb_frames,length(h.speckelQualityRoiCenterX));

                %-- Setting axis limits (mm)
                x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
                z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3; 

                %-- Setting dynamic range for visualization
                vrange = [-h.dynamic_range 0];            

                %-- padding variable
                halfWidthROIx = h.speckleQualityRoiPsfTimeX * h.system_lateralResolution;
                halfHeightROIz = h.speckleQualityRoiPsfTimeZ * h.system_axialResolution;            

                %-- Loop over frames
                for f=frame_list            

                    %-- Compute dB values
                    env = double(h.image.data(:,:,f));
                    env(env<0) = 0; env(env==0) = eps;  %-- 0's are not good friends with log10(.)
                    bmode = 20*log10(env./max(env(:)));                        
                    
                    %-- Ploting image reconstruction
                    if (h.flagDisplay==1)
                        hid = figure(1); set(gca,'fontsize',16);
                        imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                        shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                        axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                        set(gca,'YDir','reverse');
                        set(gca,'fontsize',16);
                        set(hid,'position',[124 175 1142 413]);
                        title(sprintf('Speckle Quality \n %s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                        pause(0.1);
                    end

                    %-- Main loop
                    for k=1:length(h.speckelQualityRoiCenterX)

                        %-- Compute mask inside
                        x = h.speckelQualityRoiCenterX(k);
                        z = h.speckelQualityRoiCenterZ(k);

                        %-- Compute mask ROI
                        maskROI = k * ( (h.scan.x_matrix > (x-halfWidthROIx(k))) & ...
                                     (h.scan.x_matrix < (x+halfWidthROIx(k))) & ...
                                     (h.scan.z_matrix > (z-halfHeightROIz(k))) & ...
                                     (h.scan.z_matrix < (z+halfHeightROIz(k))) );

                        %-- Extract corresponding enveloppe ROI
                        [idzz,idxx] = find(maskROI==k);
                        envRoi = env(min(idzz):max(idzz),min(idxx):max(idxx));

                        %-- Downsample the block to ensure statistical independency
                        sample = envRoi(1:5:end,1:5:end);
                        sample = sample(:);         

                        %-- Empirical evaluation of the variance on the selected block
                        var_block = tools.get_rayleigh_param(sample);

                        %-- Application of the KS test against the Rayleigh pdf with 5% confidence interval
                        score = 1-kstest(sample, 'CDF', [sample, raylcdf(sample, sqrt(var_block))], 'alpha', 0.05);

                        %-- Ploting Roi contour along with a code color
                        %-- green: test passed || red: test failed                    
                        if (h.flagDisplay==1)
                            figure(1);
                            if (score==0)
                                hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI,[1 1],'r-','linewidth',2);
                            else
                                hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI,[1 1],'g-','linewidth',2);
                            end
                            pause(0.5);
                        end
                        h.scoreSpeckleQuality(f,k) = score;                    
                    end
                    pause(2);
                end
            else
                disp('No speckle quality test')
            end
        end
        


        %------------------------------------------------------
        %-- Assess geometrical distortion
        function evaluateGeometricalDistortion(h)
  
            if (~isempty(h.geometricalDistortionX))
            
                disp('Perform geometrical distortion test')
                
                %-- Define parameters / variables
                nb_frames = length(h.image.number_plane_waves);
                frame_list = 1:nb_frames;
                h.scoreGeometricalDistortion = zeros(nb_frames,length(h.geometricalDistortionX));
                maskInside = zeros(size(h.scan.x_matrix));
                maskROI = zeros(size(h.scan.x_matrix));

                %-- Apply update on h.pht.sca -> z in order to take into
                %-- account the elevation focus
                h.geometricalDistortionZ = h.geometricalDistortionZ + h.system_z_correction;

                for k=1:length(h.geometricalDistortionX)

                    %-- Compute mask inside
                    x = h.geometricalDistortionX(k);
                    z = h.geometricalDistortionZ(k);
                    mask = k * ( (h.scan.x_matrix > (x-h.geometricalDistortionHalfWidthROI1)) & ...
                    (h.scan.x_matrix < (x+h.geometricalDistortionHalfWidthROI1)) & ...
                    (h.scan.z_matrix > (z-h.geometricalDistortionHalfWidthROI1)) & ...
                    (h.scan.z_matrix < (z+h.geometricalDistortionHalfWidthROI1)) );
                    maskInside = maskInside + mask;

                    %-- Compute mask outside
                    mask = k * ( (h.scan.x_matrix > (x-h.geometricalDistortionHalfWidthROI2)) & ...
                    (h.scan.x_matrix < (x+h.geometricalDistortionHalfWidthROI2)) & ...
                    (h.scan.z_matrix > (z-h.geometricalDistortionHalfWidthROI2)) & ...
                    (h.scan.z_matrix < (z+h.geometricalDistortionHalfWidthROI2)) );
                    maskROI = maskROI + mask;

                end

                %-- Setting axis limits (mm)
                x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
                z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3;             

                %-- Setting dynamic range for visualization
                vrange = [-h.dynamic_range 0];            

                %-- Loop over frames
                for f=frame_list            

                    %-- Compute dB values
                    env = double(h.image.data(:,:,f));
                    env(env<0) = 0; env(env==0) = eps;  %-- 0's are not good friends with log10(.)
                    bmode = 20*log10(env./max(env(:)));                        
                    
                    %-- Ploting image reconstruction
                    if (h.flagDisplay==1)
                        hid = figure(1); subplot(1,2,1); set(gca,'fontsize',16);
                        imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                        shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                        axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                        set(gca,'YDir','reverse');
                        set(gca,'fontsize',16);                           
                        title(sprintf('Geometrical Distortion \n %s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                        set(hid,'position',[124 175 1142 413]);
                    end

                    %-- Perform distortion measurements
                    for k=1:length(h.geometricalDistortionX)

                        %-- Concentrate on region of interest
                        imTest = bmode;
                        imTest(maskROI~=k) = min(bmode(:));
                        maskTestROI = maskROI;
                        maskTestROI(maskROI~=k) = 0;
                        maskTestIn = maskInside;
                        maskTestIn(maskInside~=k) = 0;

                        %-- Perform test
                        [idz,idx] = find(imTest==max(imTest(:)));    
                        if ( maskTestIn(idz,idx) == k )
                            h.scoreGeometricalDistortion(f,k) = 1;
                        end    

                        %-- Display intermediate testing
                        if (h.flagDisplay==1)                

                            figure(1); subplot(1,2,1); 
                            if (h.scoreGeometricalDistortion(f,k) == 1)
                                hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestROI,[1 1],'g-','linewidth',2);
                            else
                                hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestROI,[1 1],'r-','linewidth',2);
                            end                        
                            % hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestIn,[1 1],'g-');
                            [idzz,idxx] = find(maskTestROI==k);
                            x_lim_test = [h.scan.x_axis(min(idxx)) h.scan.x_axis(max(idxx))]*1e3;
                            z_lim_test = [h.scan.z_axis(min(idzz)) h.scan.z_axis(max(idzz))]*1e3;
                            figure(1); subplot(1,2,2); 
                            imagesc(h.scan.x_axis*1e3,h.scan.z_axis*1e3,imTest); 
                            shading flat; colormap gray; caxis(vrange); colorbar;
                            axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); 
                            set(gca,'YDir','reverse'); set(gca,'fontsize',16); axis([x_lim_test z_lim_test]);                        
                            if (h.scoreGeometricalDistortion(f,k) == 1)
                                hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestIn,[1 1],'g-','linewidth',2);   
                            else
                                hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestIn,[1 1],'r-','linewidth',2);
                            end                                             
                            hold on; plot(h.scan.x_axis(idx)*1e3,h.scan.z_axis(idz)*1e3,'ob','linewidth',2);
                            if (h.scoreGeometricalDistortion(f,k) == 1)
                                title(['Succeeded | Pt = [',num2str(round(h.scan.x_axis(idx)*1e4)/10),' , ',num2str(round(h.scan.z_axis(idz)*1e4)/10),'] mm']);
                            else
                                title(['Failed | Pt = [',num2str(round(h.scan.x_axis(idx)*1e4)/10),' , ',num2str(round(h.scan.z_axis(idz)*1e4)/10),'] mm']);
                            end    
                            pause(1);
                        end

                    end    
                    pause(1);
                end
            else
                disp('No geometrical distortion test')
            end
        end  




        %------------------------------------------------------
        %-- Assess intensity linearity
        function evaluateLinearIntensity(h)

            if (~isempty(h.linearIntensityRegionX))
            
                disp('Perform intensity linearity test')
                
                %-- Define parameters / variables
                nb_frames = length(h.image.number_plane_waves);
                frame_list = 1:nb_frames;
                h.scoreLinearIntensity = zeros(nb_frames,1);

                %-- Setting axis limits (mm)
                x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
                z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3; 

                %-- Setting the boundary of the region of interest
                minX = h.linearIntensityRegionX(1);
                maxX = h.linearIntensityRegionX(2);
                minZ = h.linearIntensityRegionZ(1);
                maxZ = h.linearIntensityRegionZ(2);   

                %-- Setting dynamic range for visualization
                vrange = [-h.dynamic_range 0];            

                %-- Loop over frames
                for f=frame_list            

                    %-- Compute dB values
                    env = double(h.image.data(:,:,f));
                    env(env<0) = 0; env(env==0) = eps;  %-- 0's are not good friends with log10(.)
                    bmode = 20*log10(env./max(env(:)));                        
                    
                    %-- Ploting image reconstruction
                    if (h.flagDisplay==1)
                        close all; hid = figure(1); subplot(1,2,1); set(gca,'fontsize',16);
                        imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                        shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                        axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                        set(gca,'YDir','reverse');
                        set(gca,'fontsize',16);
                        title(sprintf('Linear Intensity \n %s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                        set(hid,'position',[124 175 1142 413]);
                        pause(0.1);
                    end


                    %-- Compute full real y_gradient vector
                    xm = find(h.scan.x_axis>h.minX_full,1);
                    xM = find(h.scan.x_axis<h.maxX_full,1,'last');
                    zm = find(h.scan.z_axis>minZ,1);
                    zM = find(h.scan.z_axis<maxZ,1,'last');
                    y_measure_full = mean(bmode(zm:zM,xm:xM),1);
                    x_measure_full = linspace(h.minX_full,h.maxX_full,length(y_measure_full))*1e3;

                    %-- Compute theoritical y_gradient vector restricted between the area of interest, i.e. between -10 and 10 mm along x-axis
                    x_measure = linspace(minX,maxX,1000)*1e3;
                    x_ref = x_measure;
                    y_ref = linspace(0,-40,length(x_measure));
                    y_measure = interp1(x_measure_full,y_measure_full,x_measure,'spline');
                    offset = mean(y_measure-y_ref);
                    y_ref = y_ref+offset;
                    a1_ref = (y_ref(end)-y_ref(1)) / (x_measure(end)-x_measure(1));
                    a0_ref = y_ref(1)-a1_ref*x_ref(1);

                    %-- Compute the best line passing through the data (first order function y = a0 + a1x)
                    B = [ ones(size(x_measure')) , x_measure'];
                    A = (B'*B)\B'*(y_measure');
                    a0_approx = A(1);
                    a1_approx = A(2);
                    pb = [a1_approx,a0_approx];
                    y_approx = polyval(pb,x_measure);                


                    %-- Set linear intensity score
                    if (abs(a1_approx-a1_ref)<0.5)
                        h.scoreLinearIntensity(f,1) = 1;
                    end

                    if (h.flagDisplay)
                        figure(1); subplot(1,2,1);
                        hold on;
                        if (h.scoreLinearIntensity(f,1)==0)
                            plot([minX maxX]*1e3,[minZ minZ]*1e3,'r-','linewidth',2);
                            plot([minX maxX]*1e3,[maxZ maxZ]*1e3,'r-','linewidth',2);
                            plot([minX minX]*1e3,[minZ maxZ]*1e3,'r-','linewidth',2);
                            plot([maxX maxX]*1e3,[minZ maxZ]*1e3,'r-','linewidth',2);
                        else
                            plot([minX maxX]*1e3,[minZ minZ]*1e3,'g-','linewidth',2);
                            plot([minX maxX]*1e3,[maxZ maxZ]*1e3,'g-','linewidth',2);
                            plot([minX minX]*1e3,[minZ maxZ]*1e3,'g-','linewidth',2);
                            plot([maxX maxX]*1e3,[minZ maxZ]*1e3,'g-','linewidth',2);                        
                        end
                        figure(1); subplot(1,2,2);
                        plot(x_measure_full,y_measure_full,'-g','linewidth',1);
                        hold on; plot(x_measure,y_measure,'-b','linewidth',2);
                        hold on; plot(x_ref,y_ref,'-k','linewidth',2);
                        hold on; plot(x_measure,y_approx,'-r','linewidth',2);
                        set(gca,'fontsize',16);
                        if (h.scoreLinearIntensity(f,1) == 1)
                            title(sprintf('Succeeded \n Coeff ref = %1.2f \n Coeff approx = %1.2f',a1_ref,a1_approx));
                        else
                            title(sprintf('Failed \n Coeff ref = %1.2f \n Coeff approx = %1.2f',a1_ref,a1_approx));
                        end
                        hold off; pause(3);             
                    end
                end
            else
                disp('No intensity linearity test')
            end
            
        end




        %------------------------------------------------------
        %-- Assess Full Width at Half Maximum metric
        function evaluateFWHM(h)
  
            if (~isempty(h.fwhmX))

                disp('Perform fwhm measurement')
                
                %-- Define parameters / variables
                nb_frames = length(h.image.number_plane_waves);
                frame_list = 1:nb_frames;
                h.scoreFWHM = zeros(nb_frames,length(h.fwhmX));

                %-- Compute mask used to select the region of interest
                maskROI = cell(1,length(h.fwhmX));
                for k=1:length(h.fwhmX)
                    x = h.fwhmX(k);
                    z = h.fwhmZ(k);
                    maskROI{k} = k * ( (h.scan.x_matrix > (x-h.fwhmHalfWidthROI)) & ...
                                 (h.scan.x_matrix < (x+h.fwhmHalfWidthROI)) & ...
                                 (h.scan.z_matrix > (z-h.fwhmHalfWidthROI)) & ...
                                 (h.scan.z_matrix < (z+h.fwhmHalfWidthROI)) );
                end

                %-- Setting axis limits (mm)
                x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
                z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3; 

                %-- Setting dynamic range for visualization
                vrange = [-h.dynamic_range 0];

                %-- Loop over frames
                for f=frame_list

                    %-- Compute dB values
                    env = double(h.image.data(:,:,f));
                    env(env<0) = 0; env(env==0) = eps;  %-- 0's are not good friends with log10(.)
                    bmode = 20*log10(env./max(env(:)));                    
                    

                    %-- Ploting image reconstruction
                    if (h.flagDisplay==1)
                        close all;
                        hid = figure(1); subplot(2,2,1); set(gca,'fontsize',16);
                        imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                        shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                        axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                        set(gca,'YDir','reverse');
                        set(gca,'fontsize',16);                
                        title(sprintf('FWHM \n%s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                        set(hid,'position',[124 175 869 701]);
                        for k=1:length(maskROI)
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI{k},[1 1],'b-');
                        end
                    end

                    %-- Perform  resolution measurements
                    for k=1:length(h.fwhmX)

                        %-- Concentrate on region of interest
                        patchImg = double(bmode);
                        patchImg(maskROI{k}~=k) = min(bmode(:));
                        patchMask = maskROI{k};

                        %-- Extract region of interest
                        [idzz,idxx] = find(patchMask==k);
                        x_lim_patch = [h.scan.x_axis(min(idxx)) h.scan.x_axis(max(idxx))]*1e3;
                        z_lim_patch = [h.scan.z_axis(min(idzz)) h.scan.z_axis(max(idzz))]*1e3;
                        x_patch = h.scan.x_axis(min(idxx):1:max(idxx))*1e3;
                        z_patch = h.scan.z_axis(min(idzz):1:max(idzz))*1e3;

                        %-- Extract maximum point coordinates
                        [idz,idx] = find(patchImg==max(patchImg(:)));
                        idz = idz(1); idx = idx(1);
                        signalLateral = patchImg(idz,min(idxx):max(idxx));
                        signalAxial = patchImg(min(idzz):max(idzz),idx);

                        %-- Display intermediate testing
                        if (h.flagDisplay==1)
                            figure(1); subplot(2,2,1);
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI{k},[1 1],'g-','linewidth',2);
                            %-- Center display on current patch
                            figure(1); subplot(2,2,2); 
                            imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                            shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                            axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                            set(gca,'YDir','reverse');
                            set(gca,'fontsize',16);
                            hold on; plot([x_lim_patch(1) x_lim_patch(end)],[h.scan.z_axis(idz) h.scan.z_axis(idz)]*1e3,'-b','linewidth',1);
                            hold on; plot([h.scan.x_axis(idx) h.scan.x_axis(idx)]*1e3,[z_lim_patch(1) z_lim_patch(end)],'-r','linewidth',1);
                            axis([x_lim_patch z_lim_patch]);
                        end    

                        %-- Compute FWHM both on axial and lateral directions
                        [FWHM_axial] = h.ComputeFWHM(z_patch,signalAxial,3,'-r');                    
                        [FWHM_lateral] = h.ComputeFWHM(x_patch,signalLateral,4,'-b');
                        pause(1);

                        %-- Store score
                        h.scoreFWHM(f,k,1) = FWHM_axial;
                        h.scoreFWHM(f,k,2) = FWHM_lateral;

                    end
                    pause(1);
                end
            else
                disp('No fwhm measurement')
            end
        end     



        %------------------------------------------------------
        %-- Assess resolution
        function evaluateResolution(h)
            
            if (~isempty(h.resolutionAxialPointsX))

                disp('Perform resolution measurement')
                
                %-- Define parameters / variables
                nb_frames = length(h.image.number_plane_waves);
                frame_list = 1:nb_frames;
                h.scoreResolutionAxial = zeros(nb_frames,length(h.resolutionAxialPointsX));
                h.scoreResolutionLateral = zeros(nb_frames,length(h.resolutionAxialPointsX));

                %-- Generate mask to compute axial resolution
                maskAxialResolution = h.GenerateAxialResolutionMask();

                %-- Generate mask to compute lateral resolution
                maskLateralResolution = h.GenerateLateralResolutionMask();

                %-- Setting axis limits (mm)
                x_lim = [min(h.resolutionAxialPointsX)-2e-3 max(h.resolutionAxialPointsX)+2e-3]*1e3; 
                z_lim = [min(h.resolutionAxialPointsZ)-2e-3 max(h.resolutionLateralPointsZ)+2e-3]*1e3; 

                %-- Setting dynamic range for visualization
                vrange = [-h.dynamic_range 0];            

                %-- Loop over frames
                for f=frame_list

                    %-- Compute dB values
                    env = double(h.image.data(:,:,f));
                    env(env<0) = 0; env(env==0) = eps;  %-- 0's are not good friends with log10(.)
                    bmode = 20*log10(env./max(env(:)));

                    %-- Compute axial resolution
                    h.scoreResolutionAxial(f,:) = h.ComputeAxialResolution(maskAxialResolution,bmode,x_lim,z_lim,vrange,f);

                    %-- Compute lateral resolution
                    h.scoreResolutionLateral(f,:) = h.ComputeLateralResolution(maskLateralResolution,bmode,x_lim,z_lim,vrange,f);
                    pause(1);
                end
            else
                disp('No resolution measurement')
            end
            
        end


        %------------------------------------------------------
    	%-- Assess contrast
        function evaluateContrast(h)

            if (~isempty(h.constrastOcclusionCenterX))

                disp('Perfom contrast measurement')
                
                %-- Define parameters / variables
                nb_frames = length(h.image.number_plane_waves);
                frame_list = 1:nb_frames;
                h.scoreContrast = zeros(nb_frames,length(h.constrastOcclusionRadius));

                %-- Setting axis limits (mm)
                x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
                z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3; 
                x = h.scan.x_matrix;
                z = h.scan.z_matrix;            

                %-- Setting dynamic range for visualization
                vrange = [-h.dynamic_range 0];            

                %-- Loop over frames
                for f=frame_list            

                    %-- Compute dB values
                    env = double(h.image.data(:,:,f));
                    env(env<0) = 0; env(env==0) = eps;  %-- 0's are not good friends with log10(.)
                    bmode = 20*log10(env./max(env(:)));                        
                    
                    %-- Ploting image reconstruction
                    if (h.flagDisplay==1)
                        close all;
                        hid = figure(1); subplot(1,1,1); set(gca,'fontsize',16);
                        imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                        shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                        axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                        set(gca,'YDir','reverse');
                        set(gca,'fontsize',16);
                        set(hid,'position',[124 175 1142 413]);
                        title(sprintf('Contrast \n %s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                        pause(0.1);
                    end

                    %-- Main loop
                    for k=1:length(h.constrastOcclusionRadius)

                        r = h.constrastOcclusionRadius(k);
                        rin = r - h.padding * h.system_lateralResolution;
                        rout1 = r + h.padding * h.system_lateralResolution;
                        rout2 = 1.2*sqrt(rin^2+rout1^2);
                        xc = h.constrastOcclusionCenterX(k);
                        zc = h.constrastOcclusionCenterZ(k);
                        maskOcclusion = ( ((x-xc).^2 + (z-zc).^2) <= r^2);
                        maskInside = ( ((x-xc).^2 + (z-zc).^2) <= rin^2);
                        maskOutside = ( (((x-xc).^2 + (z-zc).^2) >= rout1^2) & ...
                            (((x-xc).^2 + (z-zc).^2) <= rout2^2) );

                        inside = bmode(maskInside);
                        outside = bmode(maskOutside);

                        value = 20 * log10( abs(mean(inside)-mean(outside)) / sqrt((var(inside)+var(outside))/2) );
                        h.scoreContrast(f,k) = round(value*10) / 10;

                        %-- Ploting image reconstruction
                        if (h.flagDisplay==1)
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskOcclusion,[1 1],'y-','linewidth',2);
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskInside,[1 1],'r-','linewidth',2);
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskOutside,[1 1],'g-','linewidth',2);
                            title(sprintf('Contrast (CNR = %1.2f dB)\n %s\n %d plane waves',h.scoreContrast(f,k),char(h.image.name),h.image.number_plane_waves(f)));
                            hold off; pause(1);
                        end                    
                    end      
                    pause(1);
                end
            else
                disp('No contrast measurement')
            end                
        end     

        
        %------------------------------------------------------
    	%-- Set data information        
        function set_data_information(h,info)
            %-- Set system/probe properties
            h.system_lambda = info.system_lambda;
            h.system_axialResolution = info.system_axialResolution;
            h.system_lateralResolution = info.system_lateralResolution;
            h.system_z_correction = info.system_z_correction;
            %-- Set speckle quality information
            h.speckelQualityRoiCenterX = info.speckelQualityRoiCenterX;
            h.speckelQualityRoiCenterZ = info.speckelQualityRoiCenterZ;
            h.speckleQualityRoiPsfTimeX = info.speckleQualityRoiPsfTimeX;
            h.speckleQualityRoiPsfTimeZ = info.speckleQualityRoiPsfTimeZ;
            %-- Set geometrical distortion information
            h.geometricalDistortionHalfWidthROI1 = info.system_lambda;
            h.geometricalDistortionX = info.geometricalDistortionX;
            h.geometricalDistortionZ = info.geometricalDistortionZ;
            %-- Set intensity linearity test properties
            h.linearIntensityRegionX = info.linearIntensityRegionX;
            h.linearIntensityRegionZ = info.linearIntensityRegionZ;
            %-- Set FWHM information
            h.fwhmX = info.fwhmX;
            h.fwhmZ = info.fwhmZ;
            %-- Set contrast information
            h.constrastOcclusionRadius = info.constrastOcclusionRadius;
            h.constrastOcclusionCenterX = info.constrastOcclusionCenterX;
            h.constrastOcclusionCenterZ = info.constrastOcclusionCenterZ;
            %-- Set resolution information
            h.resolutionHalfWidthROI = info.system_lambda/1.5;
            h.padResolution = info.system_lambda*1.5;
            h.resolutionLeftToRight = info.resolutionLeftToRight;
            h.resolutionAxialPointsX = info.resolutionAxialPointsX;
            h.resolutionAxialPointsZ = info.resolutionAxialPointsZ;
            h.resolutionLateralPointsX = info.resolutionLateralPointsX;
            h.resolutionLateralPointsZ = info.resolutionLateralPointsZ;
            
        end
        
        
    end
    
    
    %-- Methods call to get attribute
    methods (Access = private)
       
        %-- Compute FWHM: width at half maximum => at -6dB
        function [res] = ComputeFWHM(h,x_axis,y_signal,num,color)   
            
            %-- Perform interpolation
            coeff = 10;
            nb_sample = length(x_axis);
            nb_interp = nb_sample * coeff;
            x_interp = linspace(x_axis(1),x_axis(end),nb_interp);
            y_interp = interp1(x_axis,y_signal,x_interp);
            
            ind = find(y_interp >= (max(y_interp)-6) );
            idx1 = min(ind);
            idx2 = max(ind);
            res = x_interp(idx2) - x_interp(idx1);            
            
            %-- Display profil
            if (h.flagDisplay==1) 
                figure(1); subplot(2,2,num);
                plot(x_interp,y_interp,color,'linewidth',2);    
                hold on; plot([x_interp(idx1) x_interp(idx1)],[-100 0],'-k','linewidth',1);
                hold on; plot([x_interp(idx2) x_interp(idx2)],[-100 0],'-k','linewidth',1);            
                hold off; ylabel('Amp [dB]');
                if (num==3)
                    title(sprintf('Axial FWHM = %02.2f [mm]',res));
                    xlabel('z [mm]');
                else
                     xlabel('x [mm]');
                    title(sprintf('Lateral FWHM = %02.2f [mm]',res));
                end
            end
            
        end
        
        
        %-- Compute axial resolution
        function score = ComputeAxialResolution(h,maskROI,bmode,x_lim,z_lim,vrange,f)

            %-- Set parameter
            score = (-1)*ones(1,length(h.resolutionAxialPointsX));
            
            %-- Ploting image reconstruction
            if (h.flagDisplay==1)
                close all;
                hid = figure(1); subplot(1,2,1); set(gca,'fontsize',16);
                imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                set(gca,'YDir','reverse');
                set(gca,'fontsize',16);                
                title(sprintf('Resolution \n%s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                set(hid,'position',[124 175 842 413]);
            end

            %-- Extract maximum point coordinates of one easy case
            minVal = min(bmode(:));
          
            %-- Perform resolution measurements            
            for k=1:length(h.resolutionAxialPointsX)

                %-- Extract maximum point coordinates along x for the first mask
                patchImg = bmode;
                patchImg(maskROI{1,k}.data~=k) = minVal;
                [idz,idx] = find(patchImg==max(patchImg(:)));
                idxm = idx(1);
                idzm = idz(1);
                
                %-- Extract maximum point coordinates along x for the second mask
                patchImg = bmode;
                patchImg(maskROI{2,k}.data~=k) = minVal;
                [idz,idx] = find(patchImg==max(patchImg(:)));
                idxM = idx(1);
                idzM = idz(1);
                
                %-- Extract signal along (ptm,ptM) with interpolation          
                xm = h.scan.x_axis(idxm);
                xM = h.scan.x_axis(idxM);
                zm = h.scan.z_axis(idzm);
                zM = h.scan.z_axis(idzM);
                coeff(1) = xM-xm;
                coeff(2) = zM-zm;
                norm = sqrt(sum(coeff.^2));
                coeff = coeff ./ norm; 
                xm = xm-h.padResolution*coeff(1);
                xM = xM+h.padResolution*coeff(1);
                zm = zm-h.padResolution*coeff(2);
                zM = zM+h.padResolution*coeff(2);
                idzm = find((h.scan.z_axis>=zm),1);
                idzM = find((h.scan.z_axis<=zM),1,'last');
                z = double(h.scan.z_axis(idzm:idzM));
                coeff = 10;
                nb_interp = coeff*length(z);
                x_axis = linspace(xm,xM,nb_interp);
                z_axis = linspace(zm,zM,nb_interp);
                signal = interp2(h.scan.x_matrix,h.scan.z_matrix,bmode,x_axis,z_axis,'spline');
                
                %-- Extract the two highest local peaks
                [pks,locs] = findpeaks(signal);
                
                %-- There should be at least two peaks
                if (length(pks)>1)
                
                    [pks,indices] = sort(pks,'descend');
                    locs = locs(indices);
                    pks = pks(1:2);
                    locs = locs(1:2);

                    %-- Check whether the two highest local peaks are seperated by a point whose intensity if -3dB lower
                    if (locs(1)>locs(2))
                        x_axis_process = x_axis(locs(2):locs(1));
                        z_axis_process = z_axis(locs(2):locs(1));
                        signal_process = signal(locs(2):locs(1));
                    else
                        x_axis_process = x_axis(locs(1):locs(2));
                        z_axis_process = z_axis(locs(1):locs(2));
                        signal_process = signal(locs(1):locs(2));
                    end

                    %-- Compute the variation of intensities between the
                    %-- second maximum and the rest of the signal defined
                    %-- between the two maxima
                    dist_dB = signal(locs(2))-signal_process;
                    [diff_dB,id] = max(dist_dB);
                    ptsDispZ = z_axis_process(id);
                    ptsDispI = signal_process(id);    
                    %-- if there is a point between the two maxima  whose
                    %-- intensity is lower than 3dB => score = 1                    
                    if (~isempty(find((dist_dB>3),1)))
                        score(k) = 1e3*sqrt((x_axis_process(1)-x_axis_process(end)).^2+(z_axis_process(1)-z_axis_process(end)).^2);
                    end
                else
                    z_axis_process = z_axis;
                    signal_process = signal;
                end

                %-- Display intermediate testing
                if (h.flagDisplay==1)
                    figure(1); subplot(1,2,1);
                    mask1 = maskROI{1,k}.data;
                    mask2 = maskROI{2,k}.data;
                    hold on; plot(x_axis*1e3,z_axis*1e3,'.y')
                    if (score(k)>0)
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask1,[1 1],'g-','linewidth',2);
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask2,[1 1],'g-','linewidth',2);
                        figure(1); subplot(1,2,2);
                        plot(z_axis*1e3,signal,'b','linewidth',2);    
                        hold on; plot(z_axis(locs)*1e3,pks,'or','linewidth',2);
                        hold on; plot(z_axis_process*1e3,signal_process,'g','linewidth',2);
                        hold on; plot(ptsDispZ*1e3,ptsDispI,'or','linewidth',2);
                        title(sprintf('Axial resolution\n Passed (diff = %1.2f dB)\n Resolution = %1.2f mm',diff_dB,score(k)));
                    else
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask1,[1 1],'r-','linewidth',2);
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask2,[1 1],'r-','linewidth',2);
                        figure(1); subplot(1,2,2);
                        plot(z_axis*1e3,signal,'b','linewidth',2);    
                        hold on; plot(z_axis(locs)*1e3,pks,'or','linewidth',2);                            
                        hold on; plot(z_axis_process*1e3,signal_process,'r','linewidth',2);                        
                        title(sprintf('Axial resolution\n Failed (diff = %1.2f dB)',diff_dB));
                    end
                    hold off; ylabel('Amp [dB]');
                    xlabel('z [mm]');
                    pause(1);                   
                end
            end
        end
        
        %-- Generate mask to compute axial resolution
        function [mask] = GenerateAxialResolutionMask(h)
            
            mask = cell(2,length(h.resolutionAxialPointsX));
            if (h.resolutionLeftToRight==1)
                indice1 = 1:length(h.resolutionAxialPointsX);
                indice2 = 2:length(h.resolutionLateralPointsX);
            else
                indice1 = length(h.resolutionAxialPointsX):-1:1;
                indice2 = (length(h.resolutionLateralPointsX)-1):-1:1;
            end            
            for k=1:length(h.resolutionAxialPointsX)
                %-- Compute mask ROI
                x1 = h.resolutionAxialPointsX(indice1(k));
                %-- apply h.system_z_correction to correct the elevation focus effect for simulation data                
                z1 = h.resolutionAxialPointsZ(indice1(k))+h.system_z_correction;
                maskTmp = k * ( (h.scan.x_matrix > (x1-h.resolutionHalfWidthROI)) & ...
                             (h.scan.x_matrix < (x1+h.resolutionHalfWidthROI)) & ...
                             (h.scan.z_matrix > (z1-h.resolutionHalfWidthROI)) & ...
                             (h.scan.z_matrix < (z1+h.resolutionHalfWidthROI)) );
                mask{1,k} = struct('data',maskTmp,'x',x1,'z',z1);
                x2 = h.resolutionLateralPointsX(indice2(k));
                %-- apply h.system_z_correction to correct the elevation focus effect for simulation data
                z2 = h.resolutionLateralPointsZ(indice2(k))+h.system_z_correction;
                maskTmp = k * ( (h.scan.x_matrix > (x2-h.resolutionHalfWidthROI)) & ...
                             (h.scan.x_matrix < (x2+h.resolutionHalfWidthROI)) & ...
                             (h.scan.z_matrix > (z2-h.resolutionHalfWidthROI)) & ...
                             (h.scan.z_matrix < (z2+h.resolutionHalfWidthROI)) );
                mask{2,k} = struct('data',maskTmp,'x',x2,'z',z2);
            end
            
        end
        
        
        %-- Generate mask to compute lateral resolution
        function [mask] = GenerateLateralResolutionMask(h)
           
            mask = cell(1,length(h.resolutionLateralPointsX));
            if (h.resolutionLeftToRight==1)
                indice = length(h.resolutionLateralPointsX):-1:1;
            else
                indice = 1:1:length(h.resolutionLateralPointsX);
            end
            for k=1:length(h.resolutionLateralPointsX)
                %-- Compute mask ROI
                x = h.resolutionLateralPointsX(indice(k));
                %-- apply h.system_z_correction to correct the elevation focus effect
                z = h.resolutionLateralPointsZ(indice(k))+h.system_z_correction;
                maskTmp = k * ( (h.scan.x_matrix > (x-h.resolutionHalfWidthROI)) & ...
                             (h.scan.x_matrix < (x+h.resolutionHalfWidthROI)) & ...
                             (h.scan.z_matrix > (z-h.resolutionHalfWidthROI)) & ...
                             (h.scan.z_matrix < (z+h.resolutionHalfWidthROI)) );
                mask{1,k} = struct('data',maskTmp,'x',x,'z',z);
            end
                        
        end
        
        
        %-- Compute lateral resolution
        function score = ComputeLateralResolution(h,maskROI,bmode,x_lim,z_lim,vrange,f)

            %-- Set parameter
            score = (-1)*ones(1,length(h.resolutionLateralPointsX)-1);
            
            %-- Ploting image reconstruction
            if (h.flagDisplay==1)
                close all;
                hid = figure(1); subplot(1,2,1); set(gca,'fontsize',16);
                imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                set(gca,'YDir','reverse');
                set(gca,'fontsize',16);                
                title(sprintf('Resolution \n%s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                set(hid,'position',[124 175 842 413]);
            end

            %-- Extract maximum point coordinates of one easy case
            minVal = min(bmode(:));
          
            %-- Perform resolution measurements
            for k=1:1:(length(h.resolutionLateralPointsX)-1)

                %-- Extract maximum point coordinates along z for the first mask
                patchImg = bmode;
                patchImg(maskROI{1,k}.data~=k) = minVal;
                [idz,idx] = find(patchImg==max(patchImg(:)));
                idx1 = idx(1);
                idz1 = idz(1);
                
                %-- Extract maximum point coordinates along z for the second mask
                patchImg = bmode;
                patchImg(maskROI{1,k+1}.data~=(k+1)) = minVal;
                [idz,idx] = find(patchImg==max(patchImg(:)));
                idx2 = idx(1);
                idz2 = idz(1);
                
                %-- Extract signal along (zm,zM) with interpolation
                if (h.resolutionLeftToRight==1)
                    minX = maskROI{1,k+1}.x-h.padResolution;
                    maxX = maskROI{1,k}.x+h.padResolution;
                    idxm = idx2; idxM = idx1;
                    idzm = idz2; idzM = idz1;
                else
                    minX = maskROI{1,k}.x-h.padResolution;
                    maxX = maskROI{1,k+1}.x+h.padResolution;
                    idxm = idx1; idxM = idx2;
                    idzm = idz1; idzM = idz2;
                end
                xm = h.scan.x_axis(idxm);
                xM = h.scan.x_axis(idxM);
                zm = h.scan.z_axis(idzm);
                zM = h.scan.z_axis(idzM);
                coeff(1) = xM-xm;
                coeff(2) = zM-zm;
                norm = sqrt(sum(coeff.^2));
                if (norm==0)
                    coeff(1) = 1;
                    coeff(2) = 0;
                else
                    coeff = coeff ./ norm; 
                end
                xm = xm-h.padResolution*coeff(1);
                xM = xM+h.padResolution*coeff(1);
                zm = zm-h.padResolution*coeff(2);
                zM = zM+h.padResolution*coeff(2);                
                idxm = find((h.scan.x_axis>=xm),1);
                idxM = find((h.scan.x_axis<=xM),1,'last');
                x = double(h.scan.x_axis(idxm:idxM));
                coeff = 10;
                nb_interp = coeff*length(x);
                x_axis = linspace(xm,xM,nb_interp);
                z_axis = linspace(zm,zM,nb_interp);
                signal = interp2(h.scan.x_matrix,h.scan.z_matrix,bmode,x_axis,z_axis,'spline');
                
                %-- Extract the two highest local peaks
                [pks,locs] = findpeaks(signal);
                
                %-- There should be at least two peaks
                if (length(pks)>1)
                    
                    [pks,indices] = sort(pks,'descend');
                    locs = locs(indices);
                    pks = pks(1:2);
                    locs = locs(1:2);

                    %-- Check whether the two highest local peaks are seperated by a point whose intensity if -3dB lower                  
                    if (locs(1)>locs(2))
                        x_axis_process = x_axis(locs(2):locs(1));
                        z_axis_process = z_axis(locs(2):1:locs(1));
                        signal_process = signal(locs(2):locs(1));                        
                    else
                        x_axis_process = x_axis(locs(1):locs(2));
                        z_axis_process = z_axis(locs(1):locs(2));
                        signal_process = signal(locs(1):locs(2));                        
                    end
                    
                    %-- Compute the variation of intensities between the
                    %-- second maximum and the rest of the signal defined
                    %-- between the two maxima
                    dist_dB = signal(locs(2))-signal_process;
                    [diff_dB,id] = max(dist_dB);
                    ptsDispX = x_axis_process(id);
                    ptsDispI = signal_process(id);      
                    %-- if there is a point between the two maxima  whose
                    %-- intensity is lower than 3dB => score = 1
                    if (~isempty(find((dist_dB>3),1)))
                        score(k) = 1e3*sqrt((x_axis_process(1)-x_axis_process(end)).^2+(z_axis_process(1)-z_axis_process(end)).^2);
                    end
                else
                    x_axis_process = x_axis;
                    signal_process = signal;
                end

                %-- Display intermediate testing
                if (h.flagDisplay==1)
                    figure(1); subplot(1,2,1);
                    mask1 = maskROI{1,k}.data;
                    mask2 = maskROI{1,k+1}.data;
                    hold on; plot(x_axis*1e3,z_axis*1e3,'.y');
                    if (score(k)>0)
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask1,[1 1],'g-','linewidth',2);
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask2,[1 1],'g-','linewidth',2);
                        figure(1); subplot(1,2,2);
                        plot(x_axis*1e3,signal,'b','linewidth',2);    
                        hold on; plot(x_axis(locs)*1e3,pks,'or','linewidth',2);
                        hold on; plot(x_axis_process*1e3,signal_process,'g','linewidth',2);
                        hold on; plot(ptsDispX*1e3,ptsDispI,'or','linewidth',2);
                        title(sprintf('Lateral resolution\n Passed (diff = %1.2f dB)\n Resolution = %1.2f mm',diff_dB,score(k)));
                    else
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask1,[1 1],'r-','linewidth',2);
                        hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,mask2,[1 1],'r-','linewidth',2);
                        figure(1); subplot(1,2,2);
                        plot(x_axis*1e3,signal,'b','linewidth',2);    
                        hold on; plot(x_axis(locs)*1e3,pks,'or','linewidth',2);                            
                        hold on; plot(x_axis_process*1e3,signal_process,'r','linewidth',2);                        
                        title(sprintf('Lateral resolution\n Failed (diff = %1.2f dB)',diff_dB));
                    end
                    hold off; ylabel('Amp [dB]');
                    xlabel('x [mm]');
                    pause(1);
                end         
            end        
        end
        
    end       
    
    
    %------------------------------------------------------
    %-- Get methods

    %-- Methods call to get attribute
    methods (Access = public)
       
        function data = getNumberPlaneWavesList(h)            
            data = h.image.number_plane_waves;
        end
        
    end    
    
end


