function F = fpt_tg_fcost(p,tau,dt,T)

    dec = @(t,tau)exp(-t/tau);
    gv  = @(p,tau)(dec(p(1:end-1),tau)-dec(p(2:end),tau))/(dec(p(1),tau)-dec(p(end),tau));
    tn = length(tau);
    
    I0=0; F=0;
    for ti=1:tn
        I0 =  (2./(gv(p,tau(ti)-dt/2)+gv(p,tau(ti)+dt/2))) .* ((gv(p,tau(ti)-dt/2)-gv(p,tau(ti)+dt/2))/dt).^2;
        F  = F  + 1 ./ (tau(ti) * sqrt(sum(I0)));        
    end
    F = F/tn;



F = abs(F)^(sum(or(p<0,p>T))+1);

