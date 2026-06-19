%% Full Weak Value Analysis (Real + Imaginary)
clc; clear; close all;

%% --- CALIBRATION FOLDER FOR H AND V ---
calibration_folder = '../calibration_0907';  % Folder containing H and V files
H_pattern = '*_3_*.csv';
V_pattern = '*_48_*.csv';


%% --- DATA PARAMETERS ---
folder = './measurement5_0307_2';       
sampling_interval_ps = 4;             % Sampling interval in picoseconds
time_step = sampling_interval_ps * 1e-12;  % Convert to seconds

%% --- CALIBRATION TIME POSITIONS FOR H AND V ---
H_files = dir(fullfile(calibration_folder, H_pattern));
V_files = dir(fullfile(calibration_folder, V_pattern));

H_paths = fullfile({H_files.folder}, {H_files.name});
V_paths = fullfile({V_files.folder}, {V_files.name});

if isempty(H_paths) || isempty(V_paths)
    error('H or V calibration files not found.');
end

t_H = extract_arrival_time(H_paths, time_step);
t_V = extract_arrival_time(V_paths, time_step);


%% --- LOAD FILES ---
files     = dir(fullfile(folder, '*.csv'));
angle_map = containers.Map('KeyType','double','ValueType','any');

for k = 1:numel(files)
    fname = files(k).name;
    t     = regexp(fname, '_(\d+)_deg_', 'tokens');
    if isempty(t), continue; end
    angle = str2double(t{1}{1});
    
    if isKey(angle_map, angle)
        lst = angle_map(angle);
    else
        lst = {};
    end
    lst{end+1} = fullfile(folder, fname);
    angle_map(angle) = lst;
end

%% --- ANGLE PROCESSING ---
angles         = sort(cell2mat(keys(angle_map)))';
scaled_states  = angles * 2;
ref_idx        = find(angles == 48, 1);
if isempty(ref_idx)
    error('No 48° data point found for normalization.');
end

%% --- REAL PART: ARRIVAL TIME ANALYSIS (High-Resolution Gradient from Poly4 Fit) ---
arrival_times_ps      = nan(size(angles));
arrival_times_std_ps  = nan(size(angles));  % new: for vertical error bars

for i = 1:numel(angles)
    files_i    = angle_map(angles(i));
    all_traces = [];

    for j = 1:numel(files_i)
        opts           = detectImportOptions(files_i{j});
        opts.DataLines = [3, Inf];
        M              = readmatrix(files_i{j}, opts);
        all_traces(:,end+1) = M(:,2);
    end

    avg_trace = mean(all_traces, 2);
    time_vec  = (0:numel(avg_trace)-1) * time_step * 1e12;  % in ps

    % 1. Define 40% amplitude threshold
    threshold = 0.4 * max(avg_trace);
    idxs = find(avg_trace >= threshold);

    % 2. Define time window for poly fit
    t_window = time_vec(idxs([1 end]));
    mask     = time_vec >= t_window(1) & time_vec <= t_window(2);
    t_fit    = time_vec(mask);
    y_fit    = avg_trace(mask);

    % 3. Poly4 fit
    p = polyfit(t_fit, y_fit, 4);

    % 4. High-res time vector and evaluation
    t_hr = linspace(t_window(1), t_window(2), 10000);
    y_hr = polyval(p, t_hr);
    grad_hr = gradient(y_hr, t_hr);  % d(y)/d(t)

    % 5. Find time of max gradient
    [~, idx_max] = max(abs(grad_hr));
    arrival_single = nan(1, numel(files_i));
    for j = 1:numel(files_i)
        arrival_single(j) = extract_arrival_time({files_i{j}}, time_step);
    end

    arrival_times_ps(i)     = mean(arrival_single);
    arrival_times_std_ps(i) = std(arrival_single);

end

% Normalize to calibration
ref_arrival = arrival_times_ps(ref_idx);
norm_arrival_ps = arrival_times_ps - ref_arrival;
norm_real_weak_value = (arrival_times_ps - t_V) / (t_H - t_V);
norm_real_weak_std = arrival_times_std_ps / abs(t_H - t_V);

%norm_real_weak_value = norm_arrival_ps / max(abs(norm_arrival_ps));


%% --- PLOT REAL PART ---
figure; hold on;

% Horizontal ±0.3 degree errors
x_err = 0.3 * ones(size(scaled_states));

% Plot error bars
errorbar(scaled_states, norm_real_weak_value, ...
    norm_real_weak_std, norm_real_weak_std, ... % vertical errors
    x_err, x_err, ...                           % horizontal errors
    'o', 'Color', 'b', 'MarkerFaceColor','b', ...
    'LineWidth', 1.5, 'CapSize', 6, ...
    'DisplayName','Real Weak Value (with error)');

xlabel("\textbf{\'{E}tat d'entr\'{e}e } $|\psi(\theta)\rangle$ (d\'{e}gr\'{e}e)", ...
    'Interpreter','latex','FontSize',14);
ylabel('\textbf{Valeur faible } $\mathcal{R}$ $\langle\hat{S}_W\rangle$ (a.u.)', ...
    'Interpreter','latex','FontSize',14);
grid on;
legend('Location','best');


%% --- IMAGINARY PART: FREQUENCY ANALYSIS ---
chirp_freq_MHz = nan(size(angles));

for i = 1:numel(angles)
    files_i    = angle_map(angles(i));
    all_traces = [];

    for j = 1:numel(files_i)
        opts           = detectImportOptions(files_i{j});
        opts.DataLines = [3, Inf];
        M              = readmatrix(files_i{j}, opts);
        all_traces(:,end+1) = M(:,2);
    end

    avg_trace = mean(all_traces, 2);
    N         = numel(avg_trace);

    % FFT → power spectrum
    S = fft(avg_trace);
    P = abs(S).^2;
    f = fftshift(fftfreq(N, time_step)).';
    P = fftshift(P);
    idx = f >= 0;

    chirp_freq_MHz(i) = sum(f(idx).*P(idx)) / sum(P(idx)) * 1e-6;
end

% Normalize imaginary part to V-state (48°)
f_ref         = chirp_freq_MHz(ref_idx);
norm_freq_MHz = chirp_freq_MHz - f_ref;

%% --- COSINE FIT FOR IMAG PART ---
ft = fittype('A*(cos(B*x + C)).^2 ', ...
    'independent','x','coefficients',{'A','B','C'});
start = [ ...
    (max(norm_freq_MHz)-min(norm_freq_MHz))/2, ...
    2*pi/(max(scaled_states)-min(scaled_states)), ...
    0];
[cos_fit, gof] = fit(scaled_states, norm_freq_MHz, ft, 'StartPoint', start);

%% --- PLOT IMAGINARY PART ---
figure; hold on;
plot(scaled_states, norm_freq_MHz, 'ko', 'LineWidth',1.5, 'DisplayName','Data');
xFit = linspace(min(scaled_states), max(scaled_states), 200);
yFit = feval(cos_fit, xFit);
plot(xFit, yFit, 'r--', 'LineWidth',2, 'DisplayName','Cosine Fit');
xlabel('Input state, $|\psi(\theta)\rangle$', ...
    'Interpreter','latex','FontSize',14);
ylabel('Weak value, $\Im\langle\hat{A}_W\rangle$ (a.u.)', ...
    'Interpreter','latex','FontSize',14);
grid on;
legend('Location','best','Interpreter','none');

%% --- TABLE OUTPUT ---
T = table(angles(:), scaled_states(:), chirp_freq_MHz(:), norm_freq_MHz(:), ...
    'VariableNames',{'Angle_deg','ScaledState','Freq_MHz','NormFreq_MHz'});
disp(T);

%% --- FUNCTIONS ---

function t_arrival = extract_arrival_time(files, time_step)
    all_traces = [];

    for j = 1:numel(files)
        opts           = detectImportOptions(files{j});
        opts.DataLines = [3, Inf];
        M              = readmatrix(files{j}, opts);
        all_traces(:,end+1) = M(:,2);
    end

    avg_trace = mean(all_traces, 2);
    time_vec  = (0:numel(avg_trace)-1) * time_step * 1e12;  % in ps

    % 40% threshold and window
    threshold = 0.4 * max(avg_trace);
    idxs = find(avg_trace >= threshold);
    t_window = time_vec(idxs([1 end]));
    mask     = time_vec >= t_window(1) & time_vec <= t_window(2);
    t_fit    = time_vec(mask);
    y_fit    = avg_trace(mask);

    % Poly4 fit and high-res derivative
    p = polyfit(t_fit, y_fit, 4);
    t_hr = linspace(t_window(1), t_window(2), 10000);
    y_hr = polyval(p, t_hr);
    grad_hr = gradient(y_hr, t_hr);

    [~, idx_max] = max(abs(grad_hr));
    t_arrival = t_hr(idx_max);
end

function f = fftfreq(n, d)
    val = 1/(n*d);
    if mod(n,2)==0
        k = -n/2 : n/2-1;
    else
        k = -(n-1)/2 : (n-1)/2;
    end
    f = k * val;
end
