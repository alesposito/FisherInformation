% Photon Partitionin Theorem
% PLOS ONE
%
% A. Esposito
%
% Optimize time gates "bottom-up" iterative method
%
% Partition photons until a good F values is reached
%
%

% [p, F] = fpt_tg_bu(tau, T)
% [p, F] = fpt_tg_bu(tau, T, options)
%
% 'tau' is one value or an array of fluorescence lifetime values (e.g.,
% (0.3:0.1:3.0). Tau is used to optimize the F value numerically.
% Tau is in nanoseconds.
%
% 'T', in nanoseconds, is the period (e.g., 12.5ns for a 80MHz laser)
%
% 'p' is the partition. p provides the edges of the time gates starting
% from 0 and arriving to T. (value in nanoseconds)
%
% 'F' is the F value. If tau is an array, the F value is provided for the
% full array. Furthermore, the range over which the F value is computed is
% 25% broader in order to evaluate the performance of the system over a
% also outside tha optimization range.
% 
% 'options' is a structure providing additional information to the
% algorithm.                                                        Default
% 'options.ch_max' is the maximum channel number required           [256]
% 'options.ch_min' is the minimum channel number required           [4]
% 'options.min_df' is the minimum variation of F-value tolerated    [1e-5]
% 'options.num_it' is the maximum number of iterations              [100]
% 'options.disp'   is a logical value. If set to 1 display
%                  the plots of the optimization procedure          [1]   
%
% When a variation of at least min_df is not reached, the iterative
% algorithm stops unless che number of channels is not yet within the
% ch_max and ch_min limits. In this case, the sensitivity of the algorithm
% increases and the optomization continue further. 
% If the process in not completed in at least num_it steps, the algorithm
% stops anyway.

% EXAMPLE
%
% [p, F] = fpt_tg_bu(2, 12.5)
% [p, F] = fpt_tg_bu((.3:.1:3.0), 12.5)
%


function [p, F] = fpt_tg_bu(varargin)



    % Simple input parsing
    tau = varargin{1};
    T   = varargin{2};
    
    if nargin==2
        opt = struct;  
    elseif nargin==3
        opt = varargin{3};
    elseif nargin<2 | nargin>3    
        error('wrong number of input arguments')
    end
    
    if ~isfield(opt, 'ch_max'), opt.ch_max = 16;   end
    if ~isfield(opt, 'ch_min'), opt.ch_min = 4;    end
    if ~isfield(opt, 'min_df'), opt.min_df = 1e-5; end
    if ~isfield(opt, 'num_it'), opt.num_it = 100;  end            
    if ~isfield(opt, 'disp'),   opt.disp   = 1;    end            
    
    % Variable definitions
    dt  = 0.01; % numerical value for the computation of derivatives
    thr = 1;    % gain threshold: only channels that after splitting 
                % provide a gain grater than "thr" are actually split.
                % If no new gate provides a benefit, this values is 
                % decreased.
                % 'thr' is Fisher information

    par0 = (0:T/2:T); % Initial partition. 

    idx    = 1;
    Fcheck = [];
    Ncheck = [];
    tn     = length(tau);

    for i=1:opt.num_it

        % Evaluate the average F-value and Fisher information on the
        % initial partition
        [F0, I0] = fpt_fvalue(par0, tau, dt);
        
        % Store F and N values
        Fcheck(i) = F0;
        Ncheck(i) = size(par0,2)-1;

        % Check if to exit loop
        if Ncheck(i) >= opt.ch_max;
            display('Exiting optimization: reached maximum number of channels')
            break

        elseif i>1 & ~all(~idx)
            % this block stops the itertive process if a too low inrecement
            % of F-values is achieved. However, this is break is triggered
            % only if the partition changed (~all(~idx))
            if Fcheck(i-1)-Fcheck(i)<opt.min_df
                display('Exiting optimization: low F-value increment')
                break
            end
        end
        
        % Bild the next partition. Split all channel in two and evaluate
        % the average F-value and Fisher information on the new partition.
        par1     = sort([(par0(1:end-1)+(par0(2:end)-par0(1:end-1))/2)'; par0(1:end)'])';        
        [F1, I1] = fpt_fvalue(par1, tau, dt);
        
        % Pair the Fisher information in between the previous and new
        % partition and identify those channels that, when split, provide a
        % net inrease of Fisher information greater than 'thr'
        I1i = sum(reshape(I1,[2 size(I1,2)/2]),1);               
        idx = I1i-I0 > thr;


        if all(~idx)
            if Ncheck(i) >= opt.ch_min;
                display('Exiting optimization: minimum requested channel number reached at this sensitivity level')
                break
            else
                display('Increasing gain sensitivity: minimum requested channel numbers not reached')
                thr = thr/2;
            end 
        else
            % UPDATE PARTITION
            par0 = par1;
            par0(find(~reshape([idx; ones(size(idx))],1,[]))+1)=[];
        end
    end

    if i==opt.num_it
        warning('Reached maximum iterations. The optimization may not have converged')
    end

    figure
        subplot(1,2,1)
        plot(Ncheck)
        subplot(1,2,2)
        plot(Fcheck)

    
    % prepare outputs
    p = par0;
    F = fpt_fvalue(par0, tau, dt);