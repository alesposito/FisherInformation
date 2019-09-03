%% NOTE
%  Any type of decay can be simulated by modifying the function definition
%  block in the file fpt_fvalue.m


%% Initialization 
clear all
close all
clc

options.ch_max = 8;           % maximum number of channels
options.ch_min = 8;            % minimum number of channels
res            = 256;           % number of channels of a reference partition
tau            = (0.5:0.1:3.0); % lifetime optimization range
T              = 12.5;          % period

%% Evaluate partitions
dt         = 0.1;                                                   % lifetime variation for numerical evaluation of derivatives         
optim_par  = fpt_tg_bu(tau, T, options);                            % optimized partition
refine_par = fminsearch(@(x)fpt_tg_fcost(x,tau,dt,T),optim_par);    % refine gates
N          = length(optim_par)-1;                                   % number of optimized channels
even_par   = (0:T/N:T);                                             % even partition with N channels
dense_par  = (0:T/res:T);                                           % reference partition


%% Analyze outcomes

tn = length(tau);


% Efficiencies of partition (F^-2)
% A number between 0 and 1. 1 performs best. 0.5 mean it would require 50%
% more time to collect a high SNR as for '1'

for ti=1:tn
    FO(ti) = fpt_fvalue(optim_par,  tau(ti), dt);
    FF(ti) = fpt_fvalue(refine_par, tau(ti), dt);
    FE(ti) = fpt_fvalue(even_par,   tau(ti), dt);
    FR(ti) = fpt_fvalue(dense_par,  tau(ti), dt);
end

average_efficiency = mean([FR; FF; FO; FE]'.^-2,1);
average_efficiency = average_efficiency/average_efficiency(1);

figure
subplot(2,1,2)
bar(average_efficiency);
for i=1:4
    h=text(i,1.25,num2str(average_efficiency(i),'%2.2f'),'HorizontalAlignment','center');
end
if tn>1
    hold on
    errorbar(average_efficiency,std(([FR; FF; FO; FE]'/average_efficiency(1)).^-2,1),'xr')
end
set(gca,'Xticklabel',{[num2str(res) '-even (ref)'],[num2str(N) '-opt2'],[num2str(N) '-opt'],[num2str(N) '-even']})
title('F^-^2 - photon efficiency relative to reference partition')    
xlabel('partition')
ylabel('photon efficiency')

subplot(2,1,1)
Y0 = [0 .2 .2 0];
X  = [dense_par(1:end-1)' dense_par(1:end-1)' dense_par(2:end)' dense_par(2:end)']/T;
Y  = repmat(Y0,[size(X,1) 1]);
patch(X',Y','r');
X = [refine_par(1:end-1)' refine_par(1:end-1)' refine_par(2:end)' refine_par(2:end)']/T;
Y = repmat(Y0+0.25,[size(X,1) 1]);
patch(X',Y','r')
X = [optim_par(1:end-1)' optim_par(1:end-1)' optim_par(2:end)' optim_par(2:end)']/T;
Y = repmat(Y0+0.5,[size(X,1) 1]);
patch(X',Y','r')
X = [even_par(1:end-1)' even_par(1:end-1)' even_par(2:end)' even_par(2:end)']/T;
Y = repmat(Y0+0.75,[size(X,1) 1]);
patch(X',Y','r')
axis off
set(gca,'xlim',[-0.2 1])
text(-0.2,0.1, 'reference')
text(-0.2,0.35,'optimized-2')
text(-0.2,0.6, 'optimized')
text(-0.2,0.85, 'even')

%%
% F-value curves
tau_ref = (0.75*tau(1):tau(2)-tau(1):tau(end)*1.25);
tn = length(tau_ref);
for ti=1:tn
    FO(ti) = fpt_fvalue(optim_par,  tau_ref(ti), dt);
    FF(ti) = fpt_fvalue(refine_par, tau_ref(ti), dt);
    FE(ti) = fpt_fvalue(even_par,   tau_ref(ti), dt);
    FR(ti) = fpt_fvalue(dense_par,  tau_ref(ti), dt);
end


figure
plot(tau_ref,FO,'--r')
hold on
plot(tau_ref,FF,'r')
plot(tau_ref,FE,'b')
plot(tau_ref,FR,'k')
set(gca,'ylim',[0.9 2])
patch([tau_ref(1) tau_ref(1) tau(1) tau(1)],[0 2 2 0],[0.1 0.1 0.1],'facealpha',0.1,'linestyle','none')
patch([tau(end) tau(end) tau_ref(end) tau_ref(end) ],[0 2 2 0],[0.1 0.1 0.1],'facealpha',0.1,'linestyle','none')
xlabel('Fluorescence lifetime (ns)')
ylabel('F-value')
set(gca,'xgrid','on','ygrid','on')
legend({'optimized','optimized-2','even','reference'})