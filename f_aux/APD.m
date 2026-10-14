%% 
function [APDxx,Vpp,dVdt_max,Vr,Vp, EAD]=APD(t,v,xx)
    % APD at different repolarization % (xx, eg 90%)
    % t: vector time
    % v: vector potential (1 AP)
    % xx: repolarization
    % Other parameters: amplitude (Vpp), max upstroke (dVdt_max), Vrest, Vpeak
    % Ex.: apd90 = APD(t,v,90);  
    
    [Vp, Vp_index]=max(v); % Caution with elevated domes
    [Vr, Vr_index]=min(v(100:end));
    Vpp=Vp-Vr;
    
    % Derivative of AP
    dVdt=zeros(length(v)-1,1);
    for i=1:(length(v)-1)
        dVdt(i)=(v(i+1)-v(i))/(t(i+1)-t(i));
    end
    
    % Detect repolarization abnormalities
    t_ind_repol= find(t>100,1); 
    if sum(dVdt(t_ind_repol:end)>0.01)>1 || Vr > -70
        APDxx=NaN; Vpp=NaN; dVdt_max=NaN; Vr=NaN; Vp=NaN;
        EAD=1;
    else
        [dVdt_max, dVdt_max_index]=max(dVdt);
        t_maxderv=t(dVdt_max_index);
    
        V_xxrep=Vp-(xx/100)*Vpp;
        v_rep=v(Vp_index:end); 
        [dif,vxx_index_cut]=min(abs(v_rep-V_xxrep));
        vxx_index=vxx_index_cut+(Vp_index-1); 
        t1=t(vxx_index); v1=v(vxx_index);
        if v(vxx_index)>V_xxrep
            t2=t(vxx_index+1); v2=v(vxx_index+1);
        end
        if v(vxx_index)<=V_xxrep % solution to error when a == b
            t2=t(vxx_index-1); v2=v(vxx_index-1);
        end
        t_vxx=t1+(V_xxrep-v1)*((t2-t1)/(v2-v1));
        APDxx=t_vxx-t_maxderv;
        EAD=0;
    end
end

