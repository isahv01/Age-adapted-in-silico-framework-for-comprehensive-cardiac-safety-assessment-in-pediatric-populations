function [t_cell, StateVars_cell, currents_cell] = model_main(y0, Gx, drug, celltype, sex, ionic, phase, beats, bcl, save_beats)
% Returns cell arrays: time from 0 to bcl for each saved beat.
if nargin < 10
    save_beats = 1;
end

options = [];
keep_start = max(1, beats - save_beats + 1);
store_t = cell(save_beats, 1);
store_Y = cell(save_beats, 1);

for n = 1:beats
    [t_beat, Y_beat] = ode15s('modORd_pediatric', [0 bcl], y0, options, Gx, drug, celltype, sex, ionic, phase);
    y0 = Y_beat(end,:);
    
    if n >= keep_start
        idx = n - keep_start + 1;
        % Store relative time (starts at 0) – no offset added
        store_t{idx} = t_beat;
        store_Y{idx} = Y_beat;
    end
end

% Compute currents for each beat (same as before, but using stored time)
num_beats = save_beats;
currents_cell = cell(num_beats, 1);

for b = 1:num_beats
    t_beat = store_t{b};
    Y_beat = store_Y{b};
    npts = length(t_beat);
    
    % Preallocate current arrays
    INa = zeros(npts,1); INaL = zeros(npts,1); Ito = zeros(npts,1);
    ICaL = zeros(npts,1); IKr = zeros(npts,1); IKs = zeros(npts,1);
    IK1 = zeros(npts,1); INaCa_i = zeros(npts,1); INaCa_ss = zeros(npts,1);
    INaK = zeros(npts,1); IKb = zeros(npts,1); INab = zeros(npts,1);
    ICab = zeros(npts,1); IpCa = zeros(npts,1); JdiffCa = zeros(npts,1);
    JdiffNa = zeros(npts,1); JdiffK = zeros(npts,1); Jup = zeros(npts,1);
    Jleak = zeros(npts,1); Jtr = zeros(npts,1); Jrel = zeros(npts,1);
    CaMKa = zeros(npts,1); Istim = zeros(npts,1);
    
    for i = 1:npts
        [~, IsJs] = model_Gx_completo(t_beat(i), Y_beat(i,:), [], Gx, drug, celltype, sex, ionic, phase);
        INa(i) = IsJs(1);
        INaL(i) = IsJs(2);
        Ito(i) = IsJs(3);
        ICaL(i) = IsJs(4);
        IKr(i) = IsJs(5);
        IKs(i) = IsJs(6);
        IK1(i) = IsJs(7);
        INaCa_i(i) = IsJs(8);
        INaCa_ss(i) = IsJs(9);
        INaK(i) = IsJs(10);
        IKb(i) = IsJs(11);
        INab(i) = IsJs(12);
        ICab(i) = IsJs(13);
        IpCa(i) = IsJs(14);
        JdiffCa(i) = IsJs(15);
        JdiffNa(i) = IsJs(16);
        JdiffK(i) = IsJs(17);
        Jup(i) = IsJs(18);
        Jleak(i) = IsJs(19);
        Jtr(i) = IsJs(20);
        Jrel(i) = IsJs(21);
        CaMKa(i) = IsJs(22);
        Istim(i) = IsJs(23);
    end
    
    currents_cell{b} = [INa, INaL, Ito, ICaL, IKr, IKs, IK1, INaCa_i, INaCa_ss, ...
                         INaK, IKb, INab, ICab, IpCa, JdiffCa, JdiffNa, JdiffK, ...
                         Jup, Jleak, Jtr, Jrel, CaMKa, Istim];
end

t_cell = store_t;
StateVars_cell = store_Y;
end