function [CTDxx,cai_syst,cai_diast, EAD]=CTD(t,cai,xx, flag)  
    %Ca2+-transient duration at xx% of the initial base value (xx, eg.90%)-Passini et al.
    %For 1 APD:
    %%%% t:time vector
    %%%% cai: vector of intracellular Ca2+
    %%%% t0: stimulus time
    %%%% flag indicates if an EAD has been detected previously
    
    [cai_syst, cai_sys_index]=max(cai);
    t_ind_rep=find(t>300,1); 
    [cai_diast, cai_diast_index]=min(cai(t_ind_rep:end));
    
    dCaidt=diff(cai)./diff(t);
    
    % EADs
    if sum(dCaidt(t_ind_rep:end)>0.0001)>1 | flag == 1
        CTDxx=NaN;
        cai_syst=NaN;
        cai_diast=NaN;
        EAD = 1;
    else
        
        Acai=cai_syst-cai_diast;
        
        %CaTD
        cat_xxrec=cai_syst-(xx/100)*Acai;

        % Solution to Error: sample points must be unique!
        [cai_fragment, ia] = unique(cai(cai_sys_index:end), 'stable');
        time_i = t(cai_sys_index:end);
        time_fragment = time_i(ia);
        % Isabel HV

        t_catxx = interp1(cai_fragment,time_fragment,cat_xxrec);
        [~,ind]=max(dCaidt);
        CTDxx=t_catxx-t(ind);
        EAD = 0;
    
    end
end
