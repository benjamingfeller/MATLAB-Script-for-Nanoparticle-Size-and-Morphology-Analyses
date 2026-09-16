
function [D_A_L, D_Box] = D_finder(bw_test,edge,min_area, pixel2nm_ratio, primary_radius_nm, Rg_values_nm)
    % if edge == 1
    %     bw_test = bw_test;  %When edge = 0, imclearborder() removes any connected components touching the image border. --> as so many particles; not so important
    % elseif edge == 0          %% Problem: they are still counted in the
    %                           %%image with the red dots; so to pinpoint the number to the particle
    %                           %%this needs to be a comment
    %     bw_test = imclearborder(bw_test);
    % end
    ccBW = bwconncomp(bw_test);
    area_d = regionprops(ccBW,'Area');
    area=cell2mat(struct2cell(area_d));
    D_A_L = []; % Initialize as empty array to store valid dimensions
    D_Box = [];
    SE = strel('disk',1,0);
    for n = 1:1:ccBW.NumObjects
        if area(n) < min_area
            continue; % Skip if area is below the threshold
        end
        % Create a binary image containing only the current particle
        new_BW = zeros(ccBW.ImageSize); % Blank binary image
        new_BW(ccBW.PixelIdxList{n}) = 1; % Set current particle pixels to 1

        % Extract the bounding box around the particle/Crop image to particle bounding box (+1 pixel margin)
        [Y,X] = ind2sub([ccBW.ImageSize ccBW.ImageSize],ccBW.PixelIdxList{n});
        new_BW2 = imcrop(new_BW,[min(X)-1 min(Y)-1 max(X)-min(X)+2 max(Y)-min(Y)+2]);

        % Initialize array to store the "box area" at different resolutions
        scale_area = zeros(1, 10);
    
        % Loop through scaling factors (resolution levels)
        for i = 1:1:10;
            % Resize the particle image to simulate different box sizes using nearest-neighbour interpolation
            % example: 1 -->  11  this is upscaling factor i = 2.
            %                 11
            sBW = imresize(new_BW2,i,'Method','nearest');

            % Convert the resized image to binary (ensure consistent thresholding)
            sBW = im2bw(sBW);
            %%test:
            %figure
            %imshow(sBW)

             % Count the number of "filled boxes" (pixels) at this resolution
            scale_area(i) = length(sBW==1); % Sum of all white (1) pixels; the logical argument does not change sBW at all
        end
  
        %%eps = mean(scale_area./(1:1:10)); % Compute mean scaling factor
        %% ALL ABOVE IS NOT NEEDED!!!! IT IS A SIMPLE SOLUTION TO THE A PROP TO L^D SCALING LAW ASSUMING EVEN K=1
        eps = length(new_BW2); %%longest length of the isolated particle window/BW image
        D_A_L = [D_A_L; log(area(n)) / log(eps)]; % Append valid dimension


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % REAL BOX COUNTING (FRACTAL DIMENSION)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Convert image to logical (ensures 0/1 representation)
        % BW(y,x) = 1 → particle pixel
        % BW(y,x) = 0 → background
        BW = new_BW2 > 0;   % new_BW2 is already the croped image only containing the current particle

        % Define box sizes (ε scales)
        % Each value defines the side length of a square "box"
        %%%box_sizes = [24 28 32 36 40 44 48 52]; %--> 1 means the finest grid aka the pixels themselves; MIGHT LEAVE 1 AS TOO PIXELATED
                                  % in pixel; important: smaller boxes than
                                  % particle size! e.g. [1 4 8 12]

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OR: VARIABLE BOX SIZES (Altenhoff et al., 2020,-style BCM)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Minimum box size = average primary particle radius
        eps_min_nm = primary_radius_nm;

        % Maximum box size = radius of gyration of current aggregate
        eps_max_nm = Rg_values_nm(n);

        % Convert to pixels
        eps_min_px = eps_min_nm * pixel2nm_ratio;
        eps_max_px = eps_max_nm * pixel2nm_ratio;

        % Avoid invalid cases
        if eps_max_px <= eps_min_px
            continue
        end

        % Create 8 logarithmically spaced box sizes
        box_sizes = round(logspace( ...
            log10(eps_min_px), ...
            log10(eps_max_px), ...
            8)); %%%--> e.g. logspace(log10(10), log10(100), 8): Gives 8 values whose logarithms are evenly spaced between 1 and 2, then gives back these values in non-log space

        % Remove duplicates caused by rounding
        box_sizes = unique(box_sizes);

        % Require at least 3 scales
        if numel(box_sizes) < 3
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Will store number of occupied boxes at each scale
        Nbox = zeros(size(box_sizes));

        for k = 1:length(box_sizes)
            b = box_sizes(k);
            N = 0;

            for y = 1:b:size(BW,1)  % up to the height of the particle image
                for x = 1:b:size(BW,2)
                    block = BW( ...
                    y:min(y+b-1,end), ...   % this is a grid square! so contains multiple pixels, e.g. b=4 --> 4 x 4 pixels
                    x:min(x+b-1,end));
                    if any(block(:))  %if a block (=grid square) contains at least one 1 pixel (i.e. is (partily) part of the particle) it is counted
                        N = N + 1;    % this is the minkowski dimension
                    end
                end
            end

            Nbox(k) = N;
        end

        % log-log slope = fractal dimension
        p = polyfit(log(1./box_sizes), log(Nbox), 1);   % uses the power law N ~ eps^-D

        % fitted values
        yfit = polyval(p, log(1./box_sizes));

        % R² calculation
        SS_res = sum((log(Nbox) - yfit).^2);
        SS_tot = sum((log(Nbox) - mean(log(Nbox))).^2);
        R2 = 1 - SS_res/SS_tot;

        if n <= 5   %%% visualize first 5 box-counting fits to asses if boxes are chosen well (if a linear behavior is seen)
            figure(900 + n)
            plot(log(1./box_sizes), log(Nbox), 'bo-', 'LineWidth', 1.5); hold on;
            plot(log(1./box_sizes), yfit, 'r--', 'LineWidth', 1.5);

            xlabel('log(1/\epsilon)')
            ylabel('log(N)')
            title(sprintf('Particle %d Box Counting', n))

            grid on
            box on

            % R² annotation box
            text(min(log(1./box_sizes)), max(log(Nbox)), ...
            sprintf('R^2 = %.3f', R2), ...
            'VerticalAlignment','top', ...
            'BackgroundColor','white');

            legend('Data','Fit')
        end

        % store results (for later averaging)
        R2_all(n) = R2;
        D_all(n)  = p(1);
        
        D_Box = [D_Box; p(1)];
    end

    D_mean  = mean(D_all, 'omitnan');
    D_std   = std(D_all, 'omitnan');
    
    R2_mean = mean(R2_all, 'omitnan');
    R2_std  = std(R2_all, 'omitnan');

    %%% disp(['Mean Df Hausdorff/box counting = ' num2str(D_mean) ' ± ' num2str(D_std)]) %use the gaussian fit parameters instead of the mean
    disp(['Mean R2 as a figure of merit for box counting = ' num2str(R2_mean) ' ± ' num2str(R2_std)])
    

    %%% OPTIONAL: Cut off all to low Dfs from too small particles (when using
    %%% low minimal area threshold
    %D_Hauss = D_Hauss(D_Hauss > 1.4);



%% Plot the two different approaches for single particle Df determination:
% the first one is from the open source code, the second uses a different
% kind of scaling: the new one is a proper box counting approach: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Df from A-R exponent')
bin_edges = 1.5:0.02:max(D_A_L);
hold on
%histogram(D_A_L, bin_edges)
[counts, edges] = histcounts(D_A_L, bin_edges); %for normalised histogram
max_count = max(counts); %for normalised histogram
counts_norm = counts / max_count; %for normalised histogram
histogram('BinEdges', edges, 'BinCounts', counts_norm); %for normalised histogram

% Fit Gaussian distribution to the data
pd = fitdist(D_A_L(:), 'Normal'); % Fit normal distribution
mu = pd.mu; % Mean of the distribution
sigma = pd.sigma; % Standard deviation of the distribution
nn = numel(D_A_L); % Number of data points
error_mu = sigma / sqrt(nn); % Standard error of the mean (SEM)

% Plot Gaussian fit
x_fit = linspace(1.5, 2, 100);
y_fit = numel(D_A_L) * diff(bin_edges(1:2)) * pdf(pd, x_fit); % Scale PDF to match histogram
y_fit = y_fit / max_count; %for normalised histogram
plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Plot Gaussian fit

% Add labels, title, and legend
xlabel('Fractal Dimension (D_f)');
ylabel('Normalised Frequency');
title('Histogram D_f from A-L exponent');
lgd = legend('Histogram', 'Gaussian Fit', 'Location', 'west', 'FontSize',12);
xlim([1.55 2]);
xticks(1.55:0.05:2);
ylim([0 1.1]); %for normalised histogram
box on;
grid on;

% Display mu, sigma, and error of mu 
text(0.02, 0.98, ...
    sprintf('\\mu = %.2f, \\sigma = %.2f', mu, sigma), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',14);
hold off

figure('Name','Df from Box Counting')
bin_edges = 1.5:0.02:max(D_Box);
hold on
grid on;
box on;
%histogram(D_Box, bin_edges)
% Normalise histogram so highest bin = 1
[counts, edges] = histcounts(D_Box, bin_edges);
max_count = max(counts);

counts_norm = counts / max_count;
histogram('BinEdges', edges, 'BinCounts', counts_norm);


% Fit Gaussian distribution to the data
pd = fitdist(D_Box(:), 'Normal'); % Fit normal distribution
mu = pd.mu; % Mean of the distribution
sigma = pd.sigma; % Standard deviation of the distribution
nn = numel(D_Box); % Number of data points
error_mu = sigma / sqrt(nn); % Standard error of the mean (SEM)

% Plot Gaussian fit
x_fit = linspace(1.3, 2, 100);
y_fit = numel(D_Box) * diff(bin_edges(1:2)) * pdf(pd, x_fit); % Scale PDF to match histogram
y_fit = y_fit / max_count;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Plot Gaussian fit

% Add labels, title, and legend
xlabel('Fractal Dimension (D_f)');
ylabel('Normalised Frequency');
title('Histogram D_f from Box-Counting');
lgd = legend('Histogram', 'Gaussian Fit', 'Location', 'west', 'FontSize',12);
xlim([1.4 2]);
xticks(1.4:0.05:2);
ylim([0 1.1]);
% Set x and y font sizes.
ax = gca
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;


% Display mu, sigma, and error of mu
text(0.02, 0.98, ...
    sprintf('\\mu = %.2f, \\sigma = %.2f', mu, sigma), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',14);

hold off
disp(' ')
disp('Summary individual Df calculations using A/R and boxcounting methods:')
disp(['Mean D_log_A/R = ' num2str(mean(D_A_L))])
disp(['Mean D_Haus = ' num2str(mean(D_Box))])
disp(' ')
disp('--------------------------------------------------')



end
