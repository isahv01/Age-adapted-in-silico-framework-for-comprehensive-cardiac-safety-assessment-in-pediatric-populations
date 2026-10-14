%% 
clear, addpath('f_aux/')
%% 1. Define list of ages to make the simulations
age_list = {'infant','toddler', 'schoolAge', 'youngAdult'};

%SELECT HEART RATE:
bcl_list = {1000}; %fixed heart rate (60 bpm)
%bcl_list = {430, 500, 660, 860}; %age-adjusted heart rate

%% 2. Select options, initial conditions, drugs to simulate, and more... 
for e = 1:length(age_list)
    % load data
    age = age_list{e};
    initial_pop = load(sprintf('populations/pop_%s.mat', age));
    pop = initial_pop.poblacion(:,2:end); %from 2 because 1 is the ID
    
    % SIMULATION INPUTS
    celltype = 0;           % endo
    sex = 0;                % male
    phase = 0;              % NO hormona
    beats = 1000;           % nº of beats
    bcl = bcl_list{e};      % heart rate/cycle lentgh;
    save_beats = 3;         % nº of beats to save
    drug_list = 0;          % drug index to simulate (0-control)
    %drug_list = [1:40];    % 40 drugs from study
    d_list = 1;             % d x EFTPC -- control = 0 or 1
    %d_list = [1, 5, 10];   % (1x EFTPC, 5x and 10x) from study
    
    %Guardar simulaciones en la siguiente;
    folder = sprintf('control_sim_%s_%s', age, bcl);
    mkdir(folder)
    
    %% pararel loop
    myCluster = parcluster('local');
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj = parpool(myCluster.NumWorkers);
    end
    
    %% Select initial conditions
    if drug_list(1) == 0 %control
        y0 = getInitialVector(); %generic initial conditions
        ID_list = load(sprintf('final_300_%d', age));
    else %drug
        % if drugs are going to be simulated, initial conditions from
        % control state are needed. 
        X0 = load(sprintf('SS_final_300_%s_%d.mat', age, bcl)); SS = X0.X0;
        ID_list = SS(:, 1);
    end
    
    %% 3. BEGINNING OF THE SIMULATIONS
    % FOR EACH ID (MODEL/IONIC PROFILE)
    parfor i = 1:length(ID_list)
        ID = ID_list(i); fprintf('Simulation in progress... ID = %d\n', ID)
        
        %% FOR MAKING CONTROL SIMULATIONS KEEP THIS SECTION COMMENTED
        % % find initial contidions from control simulation
        % row_idx = find(SS(:,1) == ID, 1);
        % if isempty(row_idx)
        %     error('ID %d not found in Steady State data', ID);
        % end
        % y0 = SS(row_idx, 2:end);   % steady state variables
        %%
        % FOR EACH INDEX FROM DRUG_LIST
        for j = 1:length(drug_list)
            Num = drug_list(j);
            
            % FOR EACH CONCENTRATION
            for k = 1:length(d_list)
                d = d_list(k);
                [drug, namedrug] = drug_state(Num, d);
                profile = pop(ID, :)

                %% simulation
                [t_cell, StateVars_cell, currents_cell] = model_main(...
                    y0, profile, drug, celltype, sex, 0, phase, beats, bcl, save_beats);
                
                % create the structure to save the result
                S = struct(...
                    'ID', ID, ...
                    'pob', profile, ...
                    'namedrug', namedrug, ...
                    'd', d, ...
                    'celltype', celltype, ...
                    'sex', sex, ...
                    'phase', phase, ...
                    'edad', age, ...
                    'bcl', bcl, ...
                    't', {t_cell}, ... % to asing cell array {}
                    'StateVars', {StateVars_cell}, ...
                    'currents', {currents_cell});
                
                % SAVE SIMULATION RESULT
                mat_filename = sprintf('%s/results_ID%d_%s_%dxEFTPC_%s.mat', folder, ID, namedrug, d, age);
                parsave_struct(mat_filename, S);
            end
        end
    end
    
    % Cerrar pool paralelo
    delete(poolobj);
end

%%