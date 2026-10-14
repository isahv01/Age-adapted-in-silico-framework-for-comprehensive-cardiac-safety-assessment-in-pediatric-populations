%% 
clear, addpath('f_aux/')
%% INPUT: select from the following scenarios
edad_list = {'infant', 'toddler', 'schoolAge', 'youngAdult'};
%bcl_list = {140, 120, 90, 70}; % heart rate for each age, sorted
bcl_list = {1000, 1000, 1000, 1000}; % fixed age rate

drug_list = 1:40; % 40 drugs from study
namedrug_list = readtable('data_40_drugs.xlsx'); namedrug_list = namedrug_list.Drug;
d_list = [1, 5, 10]; % (1x EFTPC, 5x and 10x) from study

%% pararel loop
myCluster=parcluster('local');
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool(myCluster.NumWorkers);
end
%% start big loop for each age, concentration, drug and ID (model)
% FOR EACH AGE GROUP
for e = 1:length(edad_list)
    edad = edad_list{e}
    bcl = bcl_list{e};
    %
    folder = sprintf('control_sim_%s_%s', age, bcl);
    ID_list = load(sprintf('final_300_%s.csv', edad));

    temp_results = cell(25200, 12); %pre‑allocate
    final_result = []; %inicialice output
    
    %FOR EACH CONCENTRATION
    for k = 1:length(d_list)
        drug_list = d_list(k);
        
        %FOR EACH DRUG
        for i = 1:length(namedrug_list)
            namedrug = namedrug_list{i}
            
            %FOR EACH ID (MODEL) OF THE POPULATIONS
            parfor j = 1:length(ID_list)
                % select ID and load the data
                ID = ID_list(j);
                S = load(sprintf('%s/results_ID%d_%s_%dxEFTPC_%s.mat', folder, ID, namedrug, drug_list, edad))
                n_beats = length(S.t); %nº latidos guardados
                %% Compute biomarkers
                biomarks = compute_biomarkers(S, n_beats)
                
                result = {ID, namedrug, drug_list, biomarks(4), biomarks(5), biomarks(7), biomarks(12), ...
                    biomarks(18), biomarks(19), biomarks(14), biomarks(13), biomarks(18)}
                temp_results{j} = result;

            end 
            final_result = [final_result; vertcat(temp_results{:})];
        end 
    end

    %% convert to table
    T = cell2table(final_result,...
        'VariableNames', {'ID', 'DrugName', 'dxEFTPC', 'APD90_ms', ...
        'triang30_ms', 'triang50_ms', 'qnet', 'Cai_syst_µM', ...
        'Cai_diast_µM', 'CTD90_ms', 'CTD50_ms', 'EMw_ms'}); 
    % Save
    writetable(T, sprintf('biomarkers_%s_%d_w_drugs.csv', age, bcl));
end 
%% clean pararel pool
delete(poolobj); disp('Finish!!')
%