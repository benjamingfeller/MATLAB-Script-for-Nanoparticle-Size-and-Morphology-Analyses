
function [diameters, circularities] = circ_diam_finder(bw_test,edge,min_area, pixel2nm_ratio)
    if edge == 1
        bw_test = bw_test;
    elseif edge == 0
        bw_test = imclearborder(bw_test);       
    end
    
    % area determined by bwconncomp
    ccBW = bwconncomp(bw_test);
    area_d = regionprops(ccBW,'Area');
    area=cell2mat(struct2cell(area_d));

    %for better perimeter determination use same build in finder 'perim2'
    %TEST %--> Both have pros and cons;use either
    area_d2 = regionprops(ccBW,'Area', 'Perimeter');
    perim2 = [area_d2.Perimeter]'; 

    % Filter valid indices based on min_area
    valid_indices = find(area > min_area);
    area = area(valid_indices); % Filtered areas
    perim = zeros(length(valid_indices), 1); % Initialize perimeter array
    perim2 = perim2(valid_indices);  %perim with inbuilt finder now already found


    % Convert Area values to a Diameter assuming round particles:
    diameters = 2 * sqrt(area / pi);
    diameters = diameters/pixel2nm_ratio; % in nm
    % Determine circularities (use perim or perim2; see above):
    circularities = 4*pi*area(:) ./ (perim2(:).^2);  % (:) operator converts both area and perim into column vectors, avoiding dimension mismatches
    circularities = circularities';
    %perim2(:)./perim(:)

    %% Plot a histogram of diameter in log bins:

    % Define the desired bin size in log10 scale
    log_bin_size = 0.1;  % ADJUST this for finer or coarser bins

    % Compute the log10 range of diameters
    min_diameter_log = floor(log10(min(diameters)));  % Lower bound (rounded down)
    max_diameter_log = ceil(log10(max(diameters)));  % Upper bound (rounded up)

    % Create bin edges in log10 space
    bin_edges_log10 = min_diameter_log:log_bin_size:max_diameter_log;

    % Convert bin edges back to linear scale
    bin_edges = 10 .^ bin_edges_log10;

    % Plot the histogram
    figure('Name', 'Logarithmic Histogram of Diameters (OriginLab Style)');
    histogram(diameters, bin_edges, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'black');
    xlabel('Diameter (nm)');
    ylabel('Frequency');
    title(sprintf('Logarithmic Histogram of Diameters', log_bin_size)); %% (Bin Size %.1f in log10 scale)
    set(gca, 'XScale', 'log');  % Set x-axis to log scale
    grid on;

    % Define custom tick positions and labels
    custom_ticks = [0.8, 1, 2, 4, 8, 10, 20, 40, 80, 100];  % Define values for x-axis; for Au:[0.8, 1, 2, 4, 8, 10, 20, 40, 80, 100] for Pt: custom_ticks = [10, 20, 40, 80, 100];
    xticks(custom_ticks);  % Set custom tick positions
    xticklabels(string(custom_ticks));  % Set custom tick labels as their numeric values

    % Adjust the plot window x-axis range
    xlim([0.5 100]);  %for Au: 0.5-100, for Pt: xlim([8 120]);

    % Remove minor gridlines
    set(gca, 'XMinorGrid', 'off'); % Disable minor gridlines
    set(gca, 'YMinorGrid', 'off'); % Disable minor gridlines
    hold off;


    %% ================= Better alternatice: LOG-SPACE SIZE DISTRIBUTION 

% Convert diameters to log10 space
log_d = log10(diameters(:));

% Define bin width in log space
log_bin_size = 0.1;

% Calculate histogram counts
[counts, edges] = histcounts(log_d, 'BinWidth', log_bin_size);
max_count = max(counts);

% Normalise histogram so highest bin = 1
counts_norm = counts / max_count;

% Create histogram directly in log space
figure('Name','Lognormal Size Distribution');

%histogram(log_d, ...
%    'BinWidth', log_bin_size, ...
%    'FaceColor', [0.2 0.6 0.8], ...
%    'EdgeColor', 'black');

histogram('BinEdges', edges, ...
    'BinCounts', counts_norm, ...
    'FaceColor', [0.2 0.6 0.8], ...
    'EdgeColor', 'black');

xlabel('log_{10}(Diameter)');
ylabel('Normalised Frequency');
title('Particle Size Distribution');

grid on;
box on;
hold on;

%% ================= LOGNORMAL FIT (GAUSSIAN IN LOG SPACE) =================

% Fit Gaussian in log space
pd = fitdist(log_d, 'Normal');

mu_log = pd.mu;
sigma_log = pd.sigma;

n = numel(log_d);
error_mu_log = sigma_log / sqrt(n);

% Smooth x-axis in log space
x_fit = linspace(log10(1), log10(140), 300);  % adjust if another lower limit of the fit is wanted

% Gaussian fit in log space
y_fit = normpdf(x_fit, mu_log, sigma_log);

% IMPORTANT: scale to histogram counts
y_fit = y_fit * n * log_bin_size;

% Apply same normalisation as histogram
y_fit = y_fit / max_count;

% Plot fit
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

%% ================= AXIS BACK TO LINEAR VIEW (optional visual aid) =================

% original diameter scale on x-axis:
%xt = get(gca,'XTick');
%xticks(xt);
%set(gca,'XTickLabel', string(round(10.^xt)));

%%%or to get custom ticks:
xticks(log10([1 2 4 6 8 10 20 40 60])); % for Pt: xticks(log10([5 10 20 40 60 80 100 120])); Au: xticks(log10([1 2 4 6 8 10 20 40 60]));
xticklabels({'1', '2', '4', '6', '8', '10','20','40','60'}); % for Pt: xticklabels({'5','10','20','40','60','80','100','120'}); Au:xticklabels({'1', '2', '4', '6', '8', '10','20','40','60'});
xlim([log10(1) log10(60)]); % for Pt: xlim([log10(5) log10(120)]);  Au:xlim([log10(1) log10(60)]);
xlabel('Diameter (nm)');
ylim([0, 1.1]);

%% ================= ANNOTATION =================
lgd = legend('Data', 'Lognormal Fit', 'FontSize', 12);

% Display mu, sigma, and error of mu
text(0.02, 0.98, ...
    sprintf('\\mu = %.1f, \\sigma = %.1f', ...
    10^mu_log, 10^sigma_log), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',14);

hold off;

    %% Define a Logistic Model for Circularity
logistic_model = fittype('A1 - (A1 - A2) / (1 + exp((x - x0) / dx))', ...
    'independent', 'x', 'dependent', 'y');

% Perform the Fit with Initial Guesses. First make both column vectors:
diameters = diameters(:);
circularities = circularities(:);
circularities = max(circularities, 0); % Ensure all values are >= 0
circularities = min(circularities, 1); % Ensure all values are <= 1; Smallest particles overestimate the perimeter (trunkate to 1)
%%% OPTIONAL: Keep only diameters less than or equal to specific value (usefull if lots
%%% of large particles
mask = diameters <= 10;
%diameters = diameters(mask);
%circularities = circularities(mask);
%%% Make an initial guess of all parameters
initial_guess = [max(circularities), min(circularities), mean(diameters), range(diameters) / 5];
%%% where [A1, A2, x0, dx] <--> [lower asymptotic circ, higher asym ciry, midpoint, growth rate]
lower_bounds = [0, 0, mean(diameters) * 0.1, 0.1]; 
upper_bounds = [1, 1, mean(diameters) * 1.5, range(diameters)];

fit_result = fit(diameters, circularities, logistic_model, ...
    'StartPoint', initial_guess, ...
    'Lower', lower_bounds, 'Upper', upper_bounds);

%%plot and fit diameter vs circularities:
figure('Name','Circularities')
scatter(diameters,circularities, 'filled')
ylim([-inf 1.07]);
hold on;
% Plot the fitted logistic curve by evaluating the model at diameter points
x_fit = linspace(0, 20, 200);  % Fine x points for smooth curve
y_fit = feval(fit_result, x_fit);  % Evaluate the fit at the new x values
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
% Display the values of A1, A2, and x0 on the plot
A1 = fit_result.A1;
A2 = fit_result.A2;
x0 = fit_result.x0;
dx = fit_result.dx;

% Compute R^2 
y_pred = feval(fit_result, diameters);  % Predicted values from the fit
SS_res = sum((circularities - y_pred).^2);  % Residual sum of squares
SS_tot = sum((circularities - mean(circularities)).^2);  % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);
% Compute Adjusted R^2
n = length(diameters);  % Number of data points
p = 4;  % Number of fit parameters (A1, A2, x0, dx)

% Add text annotation for A1, A2, and x0, R2
str = sprintf('A1 = %.2f\nA2 = %.2f\nx_0 = %.2f\nR^2 = %.2f', ...
    A1, A2, x0, R_squared);
text(max(x_fit)-2, 0.8, str, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 12, 'BackgroundColor', 'white');

% Compute the 10% decrease point
target_circularity = A2 - 0.1 * (A2 - A1);  % Target circularity at 10% decrease
% Solve the logistic equation numerically for the corresponding diameter
logistic_eq = @(x) A1 - (A1 - A2) / (1 + exp((x - x0) / dx)) - target_circularity;

% Use fzero to find the diameter where the fit reaches target circularity
diameter_10_percent = fzero(logistic_eq, x0);  % Initial guess at midpoint x0

% Plot the point on the graph
plot(diameter_10_percent, target_circularity, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
%text(diameter_10_percent, target_circularity, sprintf('10%% decrease at %.2f nm', diameter_10_percent), ...
%    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
text(0.04, 0.98, ...
    sprintf('10%% decrease at %.2f nm', diameter_10_percent), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',13);

xlabel('Diameter (nm)');
ylabel('Circularity');
title('Logistic Fit for Circularities');
legend('Data Points', 'Logistic Fit', 'FontSize',12);
grid on;
box on;
xlim([0, 20]);
ylim([0.3, 1.1])
hold off;

%Display Fit Parameters:
%%% disp('Logistic Fit Parameters:');
%%% disp(fit_result);
%%% disp(' ');

%% plot the primary particle distribution using as a max pp size the diameter_10_percent

% Select only particles smaller than the 10% circularity decrease diameter
primary_diameters = diameters(diameters <= diameter_10_percent);

% Check that at least one particle exists
if isempty(primary_diameters)
    warning('No particles found below diameter_10_percent.');
else

    % Define the desired bin size in log10 scale
    log_bin_size = 0.1;

    % Compute log10 range
    min_diameter_log = floor(log10(min(primary_diameters)));
    max_diameter_log = ceil(log10(max(primary_diameters)));

    % Create logarithmic bin edges
    bin_edges_log10 = min_diameter_log:log_bin_size:max_diameter_log;
    bin_edges = 10.^bin_edges_log10;

    % Plot histogram
    figure('Name','Logarithmic Histogram of Primary Particle Diameters');

    histogram(primary_diameters, ...
              bin_edges, ...
              'FaceColor',[0.2 0.6 0.8], ...
              'EdgeColor','black');

    xlabel('Primary Particle Diameter (nm)');
    ylabel('Frequency');
    title('Logarithmic Histogram of Primary Particle Diameters');

    set(gca,'XScale','log');
    grid on;

    % Same tick marks as used above
    custom_ticks = [0.8, 1, 2, 4, 8, 10, 20, 40, 80, 100];
    xticks(custom_ticks);
    xticklabels(string(custom_ticks));

    xlim([0.5 100]);

    set(gca,'XMinorGrid','off');
    set(gca,'YMinorGrid','off');

end


end
