function [S_values_nm, Rg_values_nm, Df, A_values_nm, bw_img] = A_S_Rg_Dscaling_finder(gray_img, min_area, img_width_nm, particle_threshold)
    
    % Step 1: Image size and pixel-to-nm conversion
    img_width_pixels = size(gray_img, 2); % Width in pixels
    pixel2nm_ratio = img_width_nm / img_width_pixels; % Conversion factor

    % Step 2: Threshold the grayscale image to create a binary image
    % For dark particles, we need to select pixels that are less than or equal to the threshold
    bw_img = gray_img <= particle_threshold; % Pixels below the threshold are particles
    bw_img2=bw_img;
    bw_img = bwareaopen(bw_img, min_area); % Remove small particles based on area
    figure('Name', 'BW from Grayscale noise removed');
    imshow(bw_img)
    imwrite(bw_img, '225_09.png')  % if a high resolution image is needed for publications
    figure('Name', 'BW from Grayscale all pixels');
    imshow(bw_img2)

    % Step 3: Find connected components in the binary image
    ccBW = bwconncomp(bw_img);
    area_d = regionprops(ccBW, 'Area', 'PixelIdxList');
    area = [area_d.Area];

    % Step 4: Initialize arrays for S (mass-area) and Rg (radius of gyration)
    S_values = zeros(1, ccBW.NumObjects);
    Rg_values = zeros(1, ccBW.NumObjects);
    Rg_values_weighed = zeros(1, ccBW.NumObjects);
    A_values = zeros(1, ccBW.NumObjects);
    valid_idx = 0; % Counter for valid particles

    %Step 5a: Invert greyscale (such that darker means heavier/more layers of particles/thicker particles and thus has a larger value)
    inverted_img = 255 - gray_img;  % Invert pixel values (0 becomes 255, and 255 becomes 0)

    % Step 5b: Loop through connected components
    for n = 1:ccBW.NumObjects
        % Skip small particles
        if area(n) < min_area
            continue;
        end

        % Extract pixels for the current particle
        pixel_idx = ccBW.PixelIdxList{n};

        %tested but did not do a thing: Normalize the grayscale image to account for varying intensity (weighing)%
        %%normalized_img = uint8(255 * mat2gray(gray_img)); % Normalize the image intensities, better than: normalized_img = gray_img - min(gray_img(:));

        % Compute the Area (A) - number of pixels per particle
        A_values(valid_idx + 1) = length(pixel_idx);

        % Optionally subtract a background if it is significant (especially
        % if different images have differently strong background) test this
        background_intensity = mean(gray_img(:)); % Compute a global background baseline
        
        % Compute mass-area (S) as the sum of grayscale intensities
        %EITHER NO BACKGROUND SUB: 
        S_values(valid_idx + 1) = sum(inverted_img(pixel_idx)); 
        %OR WITH SUBTRACTION - test what is better(not yet a good background formula; maybe a future endeavor to continue with that but not so important):
        %S_values(valid_idx + 1) = sum(max(inverted_img(pixel_idx) - background_intensity, 0));
        % <--> Reads the pixel values of all pixels within the detected particle and sums them together.

        % Compute radius of gyration (Rg)
        [Y, X] = ind2sub(size(gray_img), pixel_idx); % Get coordinates
        cm_x = mean(X); % Center of mass x
        cm_y = mean(Y); % Center of mass y
        % Rg in pixels UNWEIGHED:
        Rg_values(valid_idx + 1) = sqrt(mean((X - cm_x).^2 + (Y - cm_y).^2)); 
        % OR Rg in pixels WEIGHED (test for S if grayscale image might be physically more relevant) normalizes intensity per pixel:
        Rg_values_weighed(valid_idx + 1) = sqrt(sum(((X - cm_x).^2 + (Y - cm_y).^2) .* double(inverted_img(pixel_idx))) / sum(double(inverted_img(pixel_idx))));

        valid_idx = valid_idx + 1;
    end

    % Keep only valid entries (not all connected components in the image
    % are valid (e.g., due to filtering by min_area). --> this steps weeds
    % them out
    S_values = S_values(1:valid_idx);
    Rg_values = Rg_values(1:valid_idx);
    Rg_values_weighed = Rg_values_weighed(1:valid_idx);
    A_values = A_values(1:valid_idx);

    % Convert Rg, A from pixels to nanometers
    Rg_values_nm = Rg_values * pixel2nm_ratio;
    Rg_values_weighed_nm = Rg_values_weighed * pixel2nm_ratio;
    A_values_nm = A_values * pixel2nm_ratio^2;
    S_values_nm = S_values * pixel2nm_ratio^2;
    

    % Step 6: Create log-log plot of S vs. Rg in nm
    log_S = log(S_values_nm);
    log_Rg = log(Rg_values_weighed_nm);  %%%OR UNWEIGHED - TEST WEIGHED HERE
    figure;
    scatter(log_Rg, log_S, 'filled');
    hold on;

    % Step 7a: Linear regression to get slope (2D fractal dimension)
    coeffs = polyfit(log_Rg, log_S, 1);
    slope = coeffs(1); % Slope of the log-log plot
    intercept = coeffs(2); % Intercept
    plot(log_Rg, polyval(coeffs, log_Rg), 'r-', 'LineWidth', 2);

    % Step 7b: Compute R^2 value and error of the slope
    y_pred = polyval(coeffs, log_Rg); % Predicted values from the fit
    SS_res = sum((log_S - y_pred).^2);  % Residual sum of squares
    SS_tot = sum((log_S - mean(log_S)).^2);  % Total sum of squares
    R_squared = 1 - (SS_res / SS_tot); % R^2

    % Strd error of the slope %% see matlab community for implementation notes
    x_var = sum((log_Rg - mean(log_Rg)).^2);  %x_var
    res = log_S - (slope*log_Rg + coeffs(2));
    res_std = std(res); %take this as an error estimate to describe the scattering of the data/scatter around the scaling law
    slope_err = res_std/sqrt(x_var*(length(log_Rg) - 1));


    % Step 8: Output results
    Df = slope; % Fractal dimension
    xlabel('log(R_g / nm)');
    ylabel('log(S / nm^2)');
    ylim([6,11]);
    yticks(6:1:11);
    xlim([-0.5,2]);
    xticks(-0.5:0.5:2);
    title('Log-Log Plot of Weighted Area S vs. R_g');
    legend('Data Points', sprintf('Fit: Slope / D_f = %.2f, R^2 = %.2f', slope, R_squared), 'Location', 'best', 'FontSize',12);
    grid on;
    box on;
    hold off;

    disp(['Bulk Fractal Dimension from Mass S vs. Radius of Gyration : ', num2str(Df), ' with error: ', num2str(res_std)]);

    % Step 9: Visualize the detected particles on the original grayscale image
    figure;
    imshow(gray_img, []); % Display the grayscale image
    hold on;

    % Loop through all detected particles and draw their boundaries
    for n = 1:ccBW.NumObjects
        if area(n) < min_area
            continue;
        end
        % Extract the pixel indices and their coordinates
        pixel_idx = ccBW.PixelIdxList{n};
        [Y, X] = ind2sub(size(gray_img), pixel_idx);

        % Compute the centroid for particle numbering
        centroid_x = mean(X);
        centroid_y = mean(Y);

        % Plot each particle's boundary in red
        plot(X, Y, 'r.', 'MarkerSize', 5); 

        % Add particle number label at the centroid
        text(centroid_x, centroid_y, num2str(n), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');

    end
    title('Detected and Analyzed Particles');
    hold off;



    %% --> Display results real A vs Rg:

% Log-Log Plot : Area (A) vs. Radius of Gyration (Rg)
log_A = log(A_values_nm);
log_Rg = log(Rg_values_nm);

figure('Name','A vs Rg')
scatter(log_Rg, log_A, 'filled');
hold on;
coeffs_A_Rg = polyfit(log_Rg, log_A, 1);
slope_A_Rg = coeffs_A_Rg(1); % Slope of the log-log plot

% Calculate the fitted values (in log space)
yfit_A_Rg = polyval(coeffs_A_Rg, log_Rg);  %function evaluated at x values log_Rg: coeffs_A_Rg(1)*log_Rg + coeffs_A_Rg(2)

% Calculate R²
SS_res_A_Rg = sum((log_A - yfit_A_Rg).^2);  % Residual sum of squares
SS_tot_A_Rg = sum((log_A - mean(log_A)).^2);  % Total sum of squares
R_squared_A_Rg = 1 - (SS_res_A_Rg / SS_tot_A_Rg);  % R² calculation

% Strd error of the slope %% see matlab community for implementation notes
x_var = sum((log_Rg - mean(log_Rg)).^2);  %x_var
res = log_A - (slope_A_Rg*log_Rg + coeffs_A_Rg(2));
res_std = std(res);  %take this as an error estimate to describe the scattering of the data/scatter around the scaling law
slope_err = res_std/sqrt(x_var*(length(log_Rg) - 1));  % only the error of the slope/ with the high number of datapoints not representative

plot(log_Rg, yfit_A_Rg, 'r-', 'LineWidth', 2);
Df_rg = slope_A_Rg; % Slope directly gives Df for area-radius scaling
xlabel('log(R_g / nm)');
ylabel('log(A / nm^2)');
title('Log-Log Plot of Area A vs. R_g ');
legend('Data Points', sprintf('Fit: Slope / D_f = %.2f, R^2 = %.2f', Df_rg, R_squared_A_Rg), 'Location', 'best', 'FontSize',12);
xlim([-0.5, 2]);
ylim([1,6]);
yticks(1:1:6);
grid on;
box on;
hold off;
disp(['Bulk Fractal Dimension from Area A vs. Radius of Gyration: ', num2str(Df_rg), ' with error: ', num2str(res_std)]);
disp(' ');

    %% test relationship Rg and Rg_weighed
    %Rg_values_nm_ratio = Rg_values_nm./Rg_values_weighed_nm
end
