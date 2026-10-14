function biomarks = compute_biomarkers(data_matfile, n_beats)
%   biomarks = compute_biomarkers(matfile, n_beats)
%
%   Inputs:
%       matfile : name of the .mat file saved by the simulation script
%       n_beats : (optional) number of beats stored (default = length(S.t))
%
%   Output:
%       biomarks : 1x18 vector of mean biomarkers (or all NaNs if any beat has EAD)
%       columns: {'APD30 (ms)','APD40 (ms)','APD50 (ms)','APD90 (ms)',...
%                   'triang30 (ms)','triang40 (ms)','triang50 (ms)','Vpp (mV)',...
%                   'dVdt_max (mV/ms)','RMP (mV)','Vpeak (mV)','qnet','CTD50 (ms)','CTD90 (ms)',...
%                   'Cai syst (µM)','Cai diast (µM)','nai (mM)', 'EMw (ms)'};
%
%   Assumes the .mat file contains fields:
%       S.t          : cell array of time vectors (one per beat)
%       S.StateVars  : cell array of state variable matrices (each row = time point)
%       S.currents   : cell array of current matrices
%
%   Requires APD() and CTD() to be in the path.

    S = data_matfile;

    % Determine number of beats from the cell array length
    if nargin < 3 || isempty(n_beats)
        n_beats = length(S.t);
    end

    % Pre‑allocate arrays
    apd30_all       = zeros(1, n_beats);
    apd40_all       = zeros(1, n_beats);
    apd50_all       = zeros(1, n_beats);
    apd90_all       = zeros(1, n_beats);
    tri30_all       = zeros(1, n_beats);
    tri40_all       = zeros(1, n_beats);
    tri50_all       = zeros(1, n_beats);
    vpp_all         = zeros(1, n_beats);        %amplitude
    dvdt_all        = zeros(1, n_beats);
    vr_all          = zeros(1, n_beats);         %RMP
    vpeak_all       = zeros(1, n_beats);
    qnet_all        = zeros(1, n_beats);     
    ctd50_all       = zeros(1, n_beats);
    ctd90_all       = zeros(1, n_beats);
    cai_syst_all    = zeros(1, n_beats);
    cai_diast_all   = zeros(1, n_beats);
    nai_all         = zeros(1, n_beats);
    EMw_all         = zeros(1, n_beats);

    % Check if currents are present
    has_currents = isfield(S, 'currents');
    
    %% start biomarkers calculation
    for i = 1:n_beats
        t = S.t{i};                         % time vector for beat i
        StateVars = S.StateVars{i};         % state variable matrix
        v = StateVars(:,1);                 % membrane voltage (mV)
        cai = StateVars(:,6) * 1e3;         % Ca²⁺ in nM (original in mM)
        nai = max(StateVars(:,2));

        % APD
        APD30 = APD(t, v, 30);
        APD40 = APD(t, v, 40);
        APD50 = APD(t, v, 50);
        [APD90, Vpp, dVdt_max, Vr, Vp, EAD] = APD(t, v, 90);
        
        tri30 = APD90 - APD30;
        tri40 = APD90 - APD40;
        tri50 = APD90 - APD50;

        % CTD (Ca²⁺ transient duration)
        [CTD90, cai_syst, cai_diast, EAD] = CTD(t, cai, 90, EAD);
        CTD50 = CTD(t, cai, 50, EAD);
        EMw = CTD90 - APD90;

        % Net charge
        currents = S.currents{i};
        % Column indices in currents matrix 
        % (assuming order from model_main.m):
        % Inet = INaL, Ito, ICaL, IKr, IKs, IK1, ... 
        Inet = currents(:,2) + currents(:,3) + ...
               currents(:,4) + currents(:,5) + ...
               currents(:,6) + currents(:,7);
        qnet = trapz(t/1000, Inet); % convert ms → s, charge in pC/pF
        qnet_all(i) = qnet;
        
        % Store results
        apd30_all(i)    = APD30;
        apd40_all(i)    = APD40;
        apd50_all(i)    = APD50;
        apd90_all(i)    = APD90;
        tri30_all(i)    = tri30;
        tri40_all(i)    = tri40;
        tri50_all(i)    = tri50;
        vpp_all(i)      = Vpp;
        dvdt_all(i)     = dVdt_max;
        vr_all(i)       = Vr;
        vpeak_all(i)    = Vp;
        ctd50_all(i)    = CTD50;
        ctd90_all(i)    = CTD90;
        cai_syst_all(i) = cai_syst;
        cai_diast_all(i)= cai_diast;
        nai_all(i)      = nai;
        EMw_all(i)      = EMw;
    end
    
    %% final vector
    % If any beat exhibited EAD, return all NaNs
    if any(EAD_flag)
        biomarks = NaN(1,10);
    else
        biomarks = [mean(apd30_all), mean(apd40_all), mean(apd50_all),...
                    mean(apd90_all), mean(tri30_all), mean(tri40_all), ...
                    mean(tri50_all), mean(vpp_all), mean(dvdt_all), ...
                    mean(vr_all), mean(vpeak_all), mean(qnet_all), ...
                    mean(ctd50_all), mean(ctd90_all), mean(cai_syst_all),...
                    mean(cai_diast_all), mean(nai_all), mean(EMw_all)];
    end
end

