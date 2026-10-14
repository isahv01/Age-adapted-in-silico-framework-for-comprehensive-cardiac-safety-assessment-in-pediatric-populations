function [t, leads] = ecg_leads_calcul(ecg_dat)
%   input = ECG (.dat from simulation)
%   output;
%   -- [t, leads] = (time, I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6)
    
    %reorganize data into variables. 
    t = ecg_dat(:, 1);
    V1 = ecg_dat(:, 2);
    V2 = ecg_dat(:, 3);
    V3 = ecg_dat(:, 4);
    V4 = ecg_dat(:, 5);
    V5 = ecg_dat(:, 6);
    V6 = ecg_dat(:, 7);
    LA = ecg_dat(:, 8);
    RA = ecg_dat(:, 9);
    LL = ecg_dat(:, 10);
    
    %calculate bipolar leads
    I = LA - RA;
    II = LL - RA;
    III = LL - LA;
    
    %calculate augmented leads 
    aVR = -(I+II/2);
    aVL = I - (II/2);
    aVF = II - (I/2);
    
    %reorganize names and variables
    leads = [I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6]; %sorted by clinical data

end
