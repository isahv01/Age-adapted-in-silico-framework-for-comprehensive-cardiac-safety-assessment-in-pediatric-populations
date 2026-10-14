clear, clc
%% CALIBRATION_0D.M  
% 1. Read the biomarkers calculated from control simulation of the complete
% population of models (10k per age).
% 2. Define the limits of the biomarkers
% 3. Do the calibration itself

%% 1. Read biomarkes data
% data_biom = {'ID', 'DrugName', 'dxEFTPC', 'APD30_ms','APD40_ms',...
    % 'APD50_ms','APD90_ms','triang30_ms','triang40_ms','triang50_ms',...
    % 'Vpp _mV','dVdt_max_mV/ms','RMP_mV','Vpeak_mV','qnet','CTD50_ms',...
    % 'CTD90_ms','Cai_syst_µM','Cai_diast_µM','nai_mM', 'EMw_ms'}

% INPUTS:    
    age = 'infant'; %select from {'infant', 'toddler', 'schoolAge', 'youngAdult'}
    bcl = 1000;
    data_biom = readmatrix(sprintf('biomarkers_%s_%d_ini_pop_control.csv', age, bcl)); 

%% 2. Define the limits of the biomarkers: from Llopis et al.2022
%APD90: ms
lim_sup_APD90 = 440;    lim_inf_APD90 = 180; 
%APD40: ms
lim_sup_APD40 = 320;    lim_inf_APD40 = 85;
%APD50: ms
lim_sup_APD50 = 350;    lim_inf_APD50 = 110;
%Tri_9040: ms
lim_sup_tri9040 = 150;  lim_inf_tri9040 = 50;
%dV/dt: mV/ms (columna 11)
lim_sup_dV = 1000;      lim_inf_dV = 150;
%Vpeak: mV (columna 13) 
lim_sup_Vp  = 55;       lim_inf_Vp = 10;
%RMP: mV (columna 12)
lim_sup_Vr = -80;       lim_inf_Vr = -95; 
%CTD50: ms (columna 15)
lim_sup_50 = 420;       lim_inf_50 = 120;
%CTD90: ms (columna 16)
lim_sup_90 = 785;       lim_inf_90 = 220;
%Cai_sys: uM (columna 17)
lim_sup_caisys = 2.23;  lim_inf_caisys = 0.262;
%Cai_dias: uM (columna 18)
lim_sup_cadias = 0.40;
%Na: mM (columna 19)
lim_sup_nai = 39.27;

%% 3. Do the calibration itself
filtered = [];
indx_filtered = [];
tam = size(data_biom);

for i = 1:tam(1)
    % Verificar todos los criterios
    cond1 = data_biom(i,7) <= lim_sup_APD90 && data_biom(i,7) >= lim_inf_APD90;     % APD90
    cond2 = data_biom(i,5) <= lim_sup_APD40 && data_biom(i,5) >= lim_inf_APD40;     % APD40
    cond3 = data_biom(i,6) <= lim_sup_APD50 && data_biom(i,6) >= lim_inf_APD50;     % APD50
    cond4 = data_biom(i,9) <= lim_sup_tri9040 && data_biom(i,9) >= lim_inf_tri9040; % triang40
    cond5 = data_biom(i,12) <= lim_sup_dV && data_biom(i,12) >= lim_inf_dV;         % dVdt_max
    cond6 = data_biom(i,14) <= lim_sup_Vp && data_biom(i,14) >= lim_inf_Vp;         % Vp
    cond7 = data_biom(i,13) <= lim_sup_Vr && data_biom(i,13) >= lim_inf_Vr;         % Vr
    cond8 = data_biom(i,16) <= lim_sup_50 && data_biom(i,16) >= lim_inf_50;         % CTD50
    cond9 = data_biom(i,17) <= lim_sup_90 && data_biom(i,17) >= lim_inf_90;         % CTD90
    cond10 = data_biom(i,18) <= lim_sup_caisys && data_biom(i,18) >= lim_inf_caisys; % cai_syst
    cond11 = data_biom(i,19) <= lim_sup_cadias;                                     % cai_diast 
    cond12 = data_biom(i,20) <= lim_sup_nai;                                        % nai 
    
    % Si cumple todos los criterios
    if cond1 && cond2 && cond3 && cond4 && cond5 && cond6 && cond7 && ...
       cond8 && cond9 && cond10 && cond11 && cond12
        filtered = [filtered; data_biom(i,:)];
        indx_filtered = [indx_filtered; i];
    else
        val = sprintf('%d, not inside ranges', i)
    end
    
end

%% 4. Save new biomarkes matrix
fprintf('=== SAVING RESULT FROM CALIBRATION_0D ===\n');
fprintf('Initial number of models: %d\n', tam(1));
fprintf('Number of models after 0D calibration: %d\n', size(filtered, 1));

writematrix(filtered, sprintf('biomarkers_%s_%d_calibration_0D_control.csv', age, bcl)) 

%% Ordenar de mayor a menor (APD90) para analizar posteriormente en 3D
%ordenar el filtro por la columna del APD90 (la 6)
%filtro = load(sprintf('N_indv_%s_datos_filtrados_1000_APDad.csv', edad)); 

%filtro_ordenado = sortrows(filtro, 6, 'ascend'); 

%writematrix(filtro_ordenado, sprintf('N_indv_%s_datos_filtrados_1000_APDad.csv', edad))

%% Quedarme con 300
data_biom = load("N_indv_toddler_datos_filtrados_1000_APDad.csv");

[idx] = find(data_biom(:, 6) >= 339 & data_biom(:, 6) <= 416);
filtered = data_biom(idx, :);

% == Quedarme con 300: RANDOM SAMPLE ==
% random = randsample(1:length(filtro), 300);
% bar = filtro(random,:);
% ID_300 = bar(:, 2)';
% filtro = data(idx, :);

% == Quedarme con 300: NORMAL SAMPLE ==
% Define the normal distribution parameters
mu = (max(filtered(:, 6)) + min(filtered(:, 6)))/2;  % center of the range
sigma = (max(filtered(:, 6)) - min(filtered(:, 6)))/6; % you choose this – e.g., (295-210)/6 ≈ 14.2 for ~99% coverage

% Compute weights: probability density of each value under N(mu,sigma)
weights = normpdf(filtered(:,6), mu, sigma);

% Sample 300 rows using weights and taking each one only once
bar = datasample(filtered, 300, 'Replace', false, 'Weights', weights);

% Extract column 2 IDs
ID_300 = bar(:, 2);
ID_300 = sort(ID_300);

%SAVE IDS
%%%%%%writematrix(ID_300, 'IDs_new_300_toddler_normal_hr60.csv')
