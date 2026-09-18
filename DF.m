
%%%%%%%%%%%%%% for every new image analysed parameters to change in the code:
% Image Name and Path
% Image width in nm
% Greyscale threshold
% Min_area threshold
% Mean TEM-derived primary particle diameter
% Number of analysed particles for primary particle size determination depending on the present amount of particles
% Size threshold of analysed particles for primary particle size determination
%%%%%%%%%%%%%%

clc, clear all, close all
addpath("INSERT PATH")
%%%read image (use grayscale here. could also insert BW. but now made as all greyscale)
gray_img = imread('IMAGE.png'); %% 
ccBW = bwconncomp(gray_img);
disp(['Image Size (width) in Pixels is: ' num2str(ccBW.ImageSize(2))]) %ccBW.ImageSize(1): Number of rows (height in pixels)., 2 = width
disp(' ')

width_image_nm = 754;  %INSERT width total image in nm (from imageJ)
pixel2nm_ratio = ccBW.ImageSize(2)/width_image_nm;  % 1 nm = this amount of pixels, (ccBW.ImageSize(2): width)
disp(['1 nm = this amount of pixels: ' num2str(pixel2nm_ratio)]) 
disp(' ')

min_area_threshold = 2^2; % INSERT Set your minimum area threshold IN NM from imagej FOR ALL BUT PP SIZE DETERMINATION
min_area_threshold = round(min_area_threshold * pixel2nm_ratio^2); %Calculates the required pixel value

%%% TEM derived mean primary diameter for variable boxing in box counting
mean_primary_diameter = 4; % in nm

%NOT useful is this at the moment as the image it uses for PP is the BW denoised:
min_area_primary_calculation = 0.67^2; %INSERT IN NM Min threshold for PP calculation - just has to be higher than the noise - can be min_area_threshld if single pp exist/are resolvable!
min_area_primary_calculation = round(min_area_primary_calculation * pixel2nm_ratio^2); %Calculates the required pixel value

num_particles = 110; % INSERT Number of particles to compute PP size - high enough to be successful - see last lines in this code

%------------also insert the greyscale image (already done above now) ----------------------
%  Input: Greyscale image, minimum area threshold, image width in nm, and grayscale threshold

%gray_img = imread('Pt-only-1L_10-2023-05-1-8bit-butsquaredORIGINAL.png'); % Replace with your grayscale image file %% Pt-only-1L_10-2023-05-1-8bit-butsquaredORIGINAL.png
min_area = min_area_threshold; % Set the minimum area threshold for valid particles
img_width_nm = width_image_nm; % Define the physical width of the image in nanometers - USE INSERT FROM ABOVE if only greyscale image is used
particle_threshold = 150; % Grayscale intensity threshold (0-255 for 8-bit images) higher = brighter pixels are also part of the particle 158 OK often
%% Show the image that will be analysed:
figure('Name', 'Original Grayscale image');
imshow(gray_img)
%% Calculate Df with greyscale image

% OPTIONAL: Call the updated function
%gray_img = 255 - gray_img; %%if BW image inserted here (only needed for 'test' image)

% calculate ensemble fractal dimensions, generate a spatial particle map showing detected particles and indices
[S_values_nm, Rg_values_nm, Df_S, A_values_nm, bw_converted] = A_S_Rg_Dscaling_finder(gray_img, min_area, img_width_nm, particle_threshold);
%%Red dots = detected particle area used in calculating S

% --> Display results real S vs Rg:
disp('number of particles detected and analysed:');
disp(length(S_values_nm));
disp('Mass-Area mean (S) value:');
disp(mean(S_values_nm));
disp('Radius of Gyration (R_g) mean value in nm:');
disp(mean(Rg_values_nm));
disp(['Fractal Dimension (D_f) from real S vs. R_g plot: ', num2str(Df_S)]);

%% from here onwards: BW image needed. --> Either from imageJ (insert on top as bw_test) OR: take concersion from grayscale (bw_converted)
%%%------Use grayscale conversion as easier produced + noise filtered out: 
bw_test = bw_converted;  %% bw_converted is already de-noised --> objects below min_area were removed in S function
%% find the area and perimeter of particles and calculate and plot the circularities and find the maximum primary particle diameter, create Particle Number-size Histogram
[diameters, circularities] = circ_diam_finder(bw_test,0, min_area_threshold, pixel2nm_ratio);


%% find box-counting fractal dimension and A-L exponent
[D_A_L, D_Box] = D_finder(bw_test,0,min_area_threshold, pixel2nm_ratio, mean_primary_diameter/2, Rg_values_nm);
%% 
%-------------------------------------------------------------------
%%%FIND PP SIZE%%%
%-------------------------------------------------------------------

%primary_sizes = compute_primary_particle_size(bw_test, min_area_primary_calculation, num_particles, pixel2nm_ratio);

