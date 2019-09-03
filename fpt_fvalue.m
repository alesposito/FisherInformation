% Evaluate the Fisher information (I) and the F-value (F) on the partition
% p, lifetimes tau using the numerical value dt to compute derivatives

function [F, I] = fpt_fvalue(p, tau, dt)

% Function definition
dec = @(t,tau)exp(-t/tau);
gv = @(p,tau)(dec(p(1:end-1),tau)-dec(p(2:end),tau))/(dec(p(1),tau)-dec(p(end),tau));
fi = @(p,tau,dt)(2./(gv(p,tau-dt/2)+gv(p,tau+dt/2))) .* ((gv(p,tau-dt/2)-gv(p,tau+dt/2))/dt).^2;
fv = @(x,I) 1 ./ (x * sqrt(sum(I)));

tn = length(tau);

I=0; F=0;
for ti=1:tn
tmp = fi(p,tau(ti),dt);
F = F + fv(tau(ti),tmp); 
I = I + tmp;
end
F = F/tn;
I = I/tn;