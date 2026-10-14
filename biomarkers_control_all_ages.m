%% 
clear, addpath('f_aux/')
%% INPUT: select from the following scenarios
%age = 'infant'; bcl = 1000;
%age = 'infant'; bcl = 430;

% age = 'toddler'; bcl = 1000;
% age = 'toddler'; bcl = 500;
% 
% age = 'schoolAge'; bcl = 1000;
% age = 'schoolAge'; bcl = 660;

age = 'youngAdult'; bcl = 1000;
%age = 'youngAdult'; bcl = 860;

n_beats = 3;

% After running the simulations
% Build a list of all result files:
result_files = dir(sprintf('sim_control_per_age/control_sim_%s_%s/results_ID*.mat', age, bcl));

% Preallocate a cell array for all results
all_results = cell(length(result_files), 21);

%% BIOMARKER CALCULATION
for k = 1:length(result_files)
    matfile = result_files(k).name;
    
    % Extract metadata from filename (optional, but you already have it inside)
    % For simplicity, we load the struct to get ID, drug, etc.
    S = load(matfile);
    
    % Compute biomarkers for this file
    biomarks = compute_biomarkers_from_beats_control(S, n_beats);
    
    % Store together with metadata
    all_results(k, :) = {S.ID, S.namedrug, S.d, biomarks(1), biomarks(2), ...
                      biomarks(3), biomarks(4), biomarks(5), biomarks(6), ...
                      biomarks(7), biomarks(8), biomarks(9), biomarks(10),...
                      biomarks(11), biomarks(12), biomarks(13), biomarks(14), ...
                      biomarks(15), biomarks(16), biomarks(17), biomarks(18)};
end
all_results = sortrows(all_results, 1);

%% Convert to table
T = cell2table(all_results(:, :), ...
    'VariableNames', {'ID', 'DrugName', 'dxEFTPC', 'APD30_ms','APD40_ms',...
    'APD50_ms','APD90_ms','triang30_ms','triang40_ms','triang50_ms',...
    'Vpp _mV','dVdt_max_mV/ms','RMP_mV','Vpeak_mV','qnet','CTD50_ms',...
    'CTD90_ms','Cai_syst_µM','Cai_diast_µM','nai_mM', 'EMw_ms'});
% Save
writetable(T, sprintf('biomarkers_%s_%d_ini_pop_control.csv', age, bcl));
