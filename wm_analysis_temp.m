% plot_interference_spectra.m
clc; clear; close all;

% ─── USER SETTINGS ────────────────────────────────────────────────────────
dataDir = './measurement5_2806';   % ← change to your FFT‐CSV folder
angles  = [3, 25, 48, 70, 93];     % the five post‑selection angles
linestyles = {'-','-','-','-','-'}; % you can vary line style if you like
colors     = lines(numel(angles)); % a nice distinct color for each

% ─── FIGURE ───────────────────────────────────────────────────────────────
figure; hold on;
for k = 1:numel(angles)
    angle = angles(k);
    % find the *_<θ>_deg_*_fft.csv file(s)
    pat = sprintf('*_%d_deg_*.csv', angle);
    D = dir(fullfile(dataDir,pat));
    if isempty(D)
        warning('No FFT file found for %d°',angle);
        continue
    end
    % if you have multiple repeats, you can average them:
    P_accum = [];
    for ii = 1:numel(D)
        M = readmatrix(fullfile(dataDir,D(ii).name));
        fMHz = M(:,1);      % assume col 1 is Hz → convert to MHz
        P    = M(:,2);           % assume col 2 is power (a.u.)
        %Fs = 1/10e-9;
        %n = length(P);
        %Y = fft(P, n);
        %f = Fs*(0:(n/2))/n;
        %P = abs(Y/sqrt(n)).^2;
        P_accum = [P_accum P];
        [power, freq] = pspectrum(P)
        %plot(f,P(1:n/2+1), 'Color', colors(k,:), ...
        plot(freq,power, 'Color', colors(k,:), ...
        'LineWidth', 2, ...
         'DisplayName', sprintf('%d°',angle));
    end
    
end

% ─── FINISH UP ────────────────────────────────────────────────────────────
%xlim([200 1000]);  % MHz, adjust to taste
%ylim([0 inf]);
xlabel('Fréquence (MHz)',       'FontSize',12);
ylabel('Spectre de puissance (a.u.)','FontSize',12);
legend('Location','northeast','FontSize',10);
grid on;
title('Interférence: Spectres vs. Angle post‑sélection','FontSize',14);
