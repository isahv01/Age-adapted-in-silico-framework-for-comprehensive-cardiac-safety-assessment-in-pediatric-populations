%% 3D calibration
% -------------------------------------
% The 3D calibration was performed using binary searches guided by the 
% APD90 limit values ​​(max and min) of the populations obtained after the
% 0D calibration.
% 
% 1. After each simulation, the QTc duration was measured and compared 
% to clinical data. 
% 2. Then, 300 models per age were selected as final populations. 

%% 1. QT CALCULATION AND CALIBRATION --------------
% here is an example
    %load ecg data
    ecg_data = load("EX_ecg_aiso_2659_infant.dat");
    
    %calculate leads from electrodes
    %ecg_leads = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', ...
                 %'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
    [ecg_time, ecg_leads] = ecg_leads_calcul(ecg_data);
    %select lead
    lead = ecg_leads(:, 1);
    
    %sample frequency (Hz)
    fs = 1000/(abs(ecg_time(1) - ecg_time(2)));
    
    %normalize and stay in baseline 0
    vn = normalize(lead,'range'); %[0 1] range
    vn = vn-vn(end); %signal baseline in 0
    lead = vn;
    
    %% find peaks (R)
    [peaks, indx_peaks] = findpeaks(lead);
    %select R_peak
    [max_peak, max_peak_indx] = max(peaks); %select maximum
    R_peak = indx_peaks(max_peak_indx);
    
    %% find peaks in a searching window from 250ms until end (T)
    [peaks_w, indx_peaks_w] = findpeaks(lead(250:end)); %select window to search for T_peak
    [max_peak_w, max_peak_w_indx] = max(peaks_w); %select maximum
    T_peak = 250 + indx_peaks_w(max_peak_w_indx); %start from window's initial value (250)
    
    %search the T final point
    T_end = visualize_TRA_method(lead, fs, T_peak, R_peak, true)
    
    %% QT and QTc calculation (ms)
    bcl = 60; %cycle length
    Q_ini = 33; %stimulus reached tissue (from simulation)
    QT = abs(Q_ini - T_end), QTc = QT*sqrt(bcl/60)
    
    figure
    grid on, hold on
    plot(ecg_time, lead, LineWidth=1.5, color = 'k')
    
    plot(T_end, lead(T_end), 'mx')
    xlim([-50 800]), xlabel('t (ms)')
    hold off
    % ---------------------------------------

% now the QTc value is compared to clinical data from Rijnbeek et al. 2001
% the new limits of APD90 after the 3D calibration are in Table 2. 

%% 2. SELECT 300 MODELS FROM THE 3D CALIBRATION 
% -------------------------------------
% To select the final 300, limits from Table 2 are needed.
%

    % INPUTS: select one
    age = 'infant'; bcl = 1000; low_limit = 351; up_limit = 417; %example
    %age = 'toddler'; bcl = 1000; low_limit = 338; up_limit = 411; 
    %age = 'schoolAge'; bcl = 1000; low_limit = 300; up_limit = 367; 
    %age = 'youngAdult'; bcl = 1000; low_limit = 270; up_limit = 350; 
    
    %%
    % read data
    data = load(sprintf('biomarkers_%s_%d_calibration_0D_control.csv', age, bcl));
    
    %select those models that meet the 3D calibration
    [idx] = find(data(:, 7) >= low_limit & data(:, 7) <= up_limit);
    filter = data(idx, :);
    
    % == Select 300: NORMAL BUT RANDOM SAMPLE ==
    % Define the normal distribution parameters
    mu = (max(filter(:, 6)) + min(filter(:, 6)))/2;     % center of the range
    sigma = (max(filter(:, 6)) - min(filter(:, 6)))/6;  % e.g., (295-210)/6 ≈ 14.2 for ~99% coverage
    
    % Compute weights: probability density of each value under N(mu,sigma)
    weights = normpdf(filter(:,6), mu, sigma);
    
    % Sample 300 rows using weights and taking each one only once
    bar = datasample(filter, 300, 'Replace', false, 'Weights', weights);
    
    % Extract column 2 IDs
    ID_300 = bar(:, 2);
    ID_300 = sort(ID_300);
    
    %% SAVE FINAL IDS
    %writematrix(ID_300, sprintf('final_300_%s.csv', age))