% SIM of one trials given parameters
%
% Luca Mazzucato March 2014

% SET OPTIONS
% ParamsRun = structure containing parameters for simulation


function [all_firings, PlotData] = fun_LIF_SIM(ParamsRun)

Theta = ParamsRun.Theta; 
Sim = ParamsRun.Sim; 
stimuli = ParamsRun.stimuli; 
%Stimulus = ParamsRun.Stimulus;
TrialStimuli = ParamsRun.TrialStimuli;
Ext = ParamsRun.Ext; 
J = ParamsRun.J; 
N_e = ParamsRun.N_e; N_i = ParamsRun.N_i; 
p = ParamsRun.p; 
He = ParamsRun.He; Hi = ParamsRun.Hi; 
tau_e = ParamsRun.tau_e; tau_i = ParamsRun.tau_i; 
tausyn_e = ParamsRun.tausyn_e; tausyn_i = ParamsRun.tausyn_i; 
tau_arp = ParamsRun.tau_arp; 
J_decision = ParamsRun.J_decision;
tau_decision = ParamsRun.tau_decision; 
tausyn_decision = ParamsRun.tausyn_decision;
N_decision = ParamsRun.N_decision; 
tau_E_r_AMPA = ParamsRun.tau_E_r_AMPA;
tau_E_d_AMPA = ParamsRun.tau_E_d_AMPA;
tau_E_r_NMDA = ParamsRun.tau_E_r_NMDA;
tau_E_d_NMDA = ParamsRun.tau_E_d_NMDA;
tau_I_r_AMPA = ParamsRun.tau_I_r_AMPA;
tau_I_d_AMPA = ParamsRun.tau_I_d_AMPA;
tau_I_r_NMDA = ParamsRun.tau_I_r_NMDA;
tau_I_d_NMDA = ParamsRun.tau_I_d_NMDA;

all_firings = [];
dt = Sim.dt_step; % time step (s)
Tseq = Sim.t_Start:dt:Sim.t_End;  

%--------------------
% PARAMETERS
%--------------------
% CELL
VEreset = He*Theta(1);
VIreset = Hi*Theta(end);
%
%----------
% STIMULUS
%----------
% build external current
% baseline current for all neurons
% add stimuli on top of baseline: for each stimulus provide 
%              - profile (perc. increase on baseline current)
%              - index of selective neurons
%
% BASELINE EXTERNAL CURRENT
mu = Ext.Mu; % mu is an (N_e+N_i)-th array
if ~isempty(stimuli)
    stim = Ext.stim;
    nstim = numel(stim); % number of stimuli in current trials
end
%----------------
% SYNAPTIC FILTER
%----------------
Tau.tausyn_e = tausyn_e; % exc synaptic time (fall time)
Tau.tausyn_i = tausyn_i; % inh synaptic time (fall time)
Tau.tausyn_decision = tausyn_decision;
Tau.tau_E_r_AMPA = tau_E_r_AMPA;
Tau.tau_E_d_AMPA = tau_E_d_AMPA;
Tau.tau_E_r_NMDA = tau_E_r_NMDA;
Tau.tau_E_d_NMDA = tau_E_d_NMDA;
Tau.tau_I_r_AMPA = tau_I_r_AMPA;
Tau.tau_I_d_AMPA = tau_I_d_AMPA;
Tau.tau_I_r_NMDA = tau_I_r_NMDA;
Tau.tau_I_d_NMDA = tau_I_d_NMDA;
F = synaptic_trace(Tau,dt,N_e,N_i,N_decision,'sexp'); % traces for recurrent connections

%--------------------------
% SIMULATION
%--------------------------
% preallocate memory for stored variable firings_tmp
% INITIAL CONDITIONS: random
v = [(Theta(1)-VEreset)/2*ones(N_e,1)+(Theta(1)-VEreset)/2*(2*rand(N_e,1)-1); ...
     (Theta(end)-VIreset)/2*ones(N_i,1)+(Theta(end)-VIreset)/2*(2*rand(N_i,1)-1); ...
      VEreset*ones(2*N_decision,1)]; % Initial values of v
% THRESHOLD VECTOR
VTh = [Theta(1)*ones(N_e,1) ; Theta(end)*ones(N_i,1) ; Theta(1)*ones(2*N_decision,1)];
c = [VEreset*ones(N_e,1) ; VIreset*ones(N_i,1) ; VEreset*ones(2*N_decision,1)]; % reset potentials
% fprintf('\nVEThres = %g; VIThres = %g',VTh(1),VTh(end));
% fprintf('\nVEreset = %g; VIreset = %g \n',c(1),c(end));
     % Excitatory neurons       Inhibitory neurons     Decision neurons
tau = [tau_e*ones(N_e,1);       tau_i*ones(N_i,1);     tau_decision*ones(2*N_decision,1)];
%
firings = zeros(10*numel(Tseq),2);
firings_cnt = 0;
tic
%--------------------
% PLOT
%--------------------
PlotData = [];
PlotData.Ne_plot = N_e; % number of exc neuron to plot
PlotData.Ni_plot = N_i; % number of inh neurons to plot
ind_plot = [5 ; N_e+5]; % indices of neurons to plot
if ~isempty(stimuli)
    indcue = find(cellfun(@(x) ~isempty(x), strfind({stim(:).name},'CS')));
    if ~isempty(indcue)
        ind_plot(1) = stim(indcue).ind(1);
    end
end
nplot = numel(ind_plot); % number of neurons to plot (membrane potential plot)
vi = 0; % running index for vplot
PlotData.vplot = zeros(nplot,round(Sim.plot_length/dt)); % store membrane potential for plots; rows=neurons, cols=time steps;
PlotData.iEplot = zeros(2,round(Sim.plot_length/dt)); % store EPSC for plots; rows=neurons, cols=time steps;
PlotData.iExtplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
PlotData.iIplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
PlotData.p = p;
PlotData.VTh = VTh;
PlotData.tau = tau;
PlotData.ind_plot = ind_plot;
%----------------------------
% RUN
%----------------------------

refr = zeros(size(mu,1)+2*N_decision, 1); % neurons in refractory state
for t = 1:numel(Tseq) % simulation of 1000 ms 
    fired = find(v > VTh); % indices of spikes
    Isyn = zeros(N_e+N_i+2*N_decision, 1);
    % spikes
    if ~isempty(fired)      
        v(fired) = c(fired);  
        refr(fired) = tau_arp;
    end
    % recurrent synaptic current
    F = syn_evolve(F,fired,'sexp');    
    % integrate
    muRun = mu;
    if ~isempty(stimuli)
        for n = 1:nstim
            if ~ismember(stim(n).name,TrialStimuli), continue; end
            if Tseq(t) >= stim(n).interval(1) && Tseq(t) <= stim(n).interval(2) 
                if strcmp(stim(n).name,'CSgauss')
                    muRun(stim(n).ind) = muRun(stim(n).ind) + stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
                else
                    if any(strcmp(stim(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'}))
                        % Some neurons are excited by stimulus
                        muRun(stim(n).ind_inc) = muRun(stim(n).ind_inc) + stim(n).profile_inc(Tseq(t))*mu(stim(n).ind_inc).*stim(n).DiffGain_Inc;
                        % and others are inhibited
                        muRun(stim(n).ind_dec) = muRun(stim(n).ind_dec) + stim(n).profile_dec(Tseq(t))*mu(stim(n).ind_dec).*stim(n).DiffGain_Dec;
                    else
                        muRun(stim(n).ind) = muRun(stim(n).ind) + stim(n).profile(Tseq(t))*mu(stim(n).ind);
                    end
                end
            end
        end
    end
    if isstruct(J_decision)
        switch TrialStimuli{1}
            case 'Sucrose'
                J_D = J_decision.S;
            case 'Quinine'
                J_D = J_decision.Q;
            case 'Maltose'
                J_D = J_decision.M;
            case 'Octaacetate'
                J_D = J_decision.O;
            otherwise
                J_D = zeros(size(J_decision.S));    
        end
    else
        J_D = J_decision;
    end
    Isyn(1:(N_e+N_i)) = Isyn(1:(N_e+N_i)) + J*F.f(1:(N_e+N_i));
    Isyn((N_e+N_i+1):(N_e+N_i+2*N_decision)) = Isyn((N_e+N_i+1):(N_e+N_i+2*N_decision)) + J_D*F.f;
    v = v - v*dt./tau + [muRun(:) ; zeros(2*N_decision,1)]*dt + Isyn*dt;
    % neurons in refractory state
    refr = max(-0.001,refr-dt); 
    v(refr>0) = c(refr>0);
    % store spikes
    if ~isempty(fired)
        % if firings_tmp has no more space, preallocate more memory
        if firings_cnt + numel(fired) > size(firings,1)
            firings = [firings; zeros(10*numel(Tseq),2)];
        end
        firings(firings_cnt+1:firings_cnt+numel(fired),1:2) = [Tseq(t)+0*fired, fired];  
        firings_cnt = firings_cnt + numel(fired);
    end
    % store values for plotting, only last Sim.plot_length interval
    if Tseq(t) > Sim.t_End - Sim.plot_length
        vi = vi + 1;
        % membrane potential
        PlotData.vplot(1:nplot,vi) = v(ind_plot); 
        % raw PSC output (other neurons see a weighted sum of this and other inputs)
        PlotData.PSCplot(1:nplot,vi) = F.f(ind_plot);
        PlotData.PSCplot_dec(1:2,vi) = F.f([N_e+N_i+5 ; N_e+N_i+N_decision+5]);
        if isfield(F,'f_AMPA')
            PlotData.PSCplot_AMPA(1:nplot,vi) = F.f_AMPA(ind_plot);
            PlotData.PSCplot_NMDA(1:nplot,vi) = F.f_NMDA(ind_plot);
        end
        % input currents
        PlotData.iEplot(1:nplot,vi) = J(ind_plot,1:N_e)*F.f(1:N_e);
        PlotData.iIplot(1:nplot,vi) = J(ind_plot,N_e+1:N_e+N_i)*F.f(N_e+1:N_e+N_i);
        PlotData.iExtplot(1:nplot,vi) = muRun(ind_plot,1);
        PlotData.iE2Dplot(1:2,vi) = J_D([5 ; N_decision+5],1:N_e+N_i)*F.f(1:N_e+N_i);
        PlotData.iD2Dplot(1:2,vi) = J_D([5 ; N_decision+5],N_e+N_i+1:N_e+N_i+2*N_decision)*F.f((N_e+N_i+1):(N_e+N_i+2*N_decision));
    end
end

% fprintf('--- End of trial...\n');
% toc
%---------------------------------------------------------------------------
if ~any(any(firings))
    fprintf('\n --- NO SPIKES GENERATED... \n');
else
    % find last spike in firings
    IndexEnd = find(firings(:,2)==0,1) - 1;
    if isempty(IndexEnd)
        IndexEnd = size(firings,1);
    end
    all_firings = firings(1:IndexEnd,[1 2]);
end

%% Sub-functions
function F = synaptic_trace(Tau,dt,N_e,N_i,N_decision,dynamicsModel)
    F = struct();
    if strcmp(dynamicsModel,'sexp')
        tau_sE = Tau.tausyn_e; % exc synaptic time (fall time)
        tau_sI = Tau.tausyn_i; % inh synaptic time (fall time)
        tau_sD = Tau.tausyn_decision;
        fexp = [repmat(exp(-dt/tau_sE),N_e,1); repmat(exp(-dt/tau_sI),N_i,1); repmat(exp(-dt/tau_sD),2*N_decision,1)]; % Multiplicative step (fp)
        fSpike = [repmat((1/tau_sE),N_e,1); repmat((1/tau_sI),N_i,1) ; repmat((1/tau_sD),2*N_decision,1)]; % add to fp with a spike
        f = zeros(N_e+N_i+2*N_decision,1);
        F.fexp = fexp;
        F.fSpike = fSpike;
        F.f = f;
    elseif strcmp(dynamicsModel,'dexp')
        tau_E_r_AMPA = Tau.tau_E_r_AMPA;
        tau_E_d_AMPA = Tau.tau_E_d_AMPA;
        tau_E_r_NMDA = Tau.tau_E_r_NMDA;
        tau_E_d_NMDA = Tau.tau_E_d_NMDA;
        tau_I_r_AMPA = Tau.tau_I_r_AMPA;
        tau_I_d_AMPA = Tau.tau_I_d_AMPA;
        tau_I_r_NMDA = Tau.tau_I_r_NMDA;
        tau_I_d_NMDA = Tau.tau_I_d_NMDA;
        fexp_r_AMPA = [repmat(exp(-dt/tau_E_r_AMPA),N_e,1); repmat(exp(-dt/tau_I_r_AMPA),N_i,1)]; % Multiplicative step (fp)
        fexp_d_AMPA = [repmat(exp(-dt/tau_E_d_AMPA),N_e,1); repmat(exp(-dt/tau_I_d_AMPA),N_i,1)];
        fexp_r_NMDA = [repmat(exp(-dt/tau_E_r_NMDA),N_e,1); repmat(exp(-dt/tau_I_r_NMDA),N_i,1)];
        fexp_d_NMDA = [repmat(exp(-dt/tau_E_d_NMDA),N_e,1); repmat(exp(-dt/tau_I_d_NMDA),N_i,1)];
        fScale_AMPA = [repmat(1/(tau_E_d_AMPA-tau_E_r_AMPA),N_e,1); repmat(1/(tau_I_d_AMPA-tau_I_r_AMPA),N_i,1)];
        fScale_NMDA = [repmat(1/(tau_E_d_NMDA-tau_E_r_NMDA),N_e,1); repmat(1/(tau_I_d_NMDA-tau_I_r_NMDA),N_i,1)];
        fSpike_r = [repmat((-1),N_e,1); repmat((-1),N_i,1)]; % add to fp with a spike
        fSpike_d = [repmat((1),N_e,1); repmat((1),N_i,1)];
        f = zeros(N_e+N_i,1);
        F.fexp_r_AMPA = fexp_r_AMPA;
        F.fexp_d_AMPA = fexp_d_AMPA;
        F.fexp_r_NMDA = fexp_r_NMDA;
        F.fexp_d_NMDA = fexp_d_NMDA;
        F.fScale_AMPA = fScale_AMPA;
        F.fScale_NMDA = fScale_NMDA;
        F.fSpike_r = fSpike_r;
        F.fSpike_d = fSpike_d;
        F.f_r_AMPA = f;
        F.f_d_AMPA = f;
        F.f_AMPA = f;
        F.f_r_NMDA = f;
        F.f_d_NMDA = f;
        F.f_NMDA = f;
        F.f = f;
    end
end

function F = syn_evolve(F,fired,dynamicsModel)
    % update synaptic filter
    if strcmp(dynamicsModel,'sexp')
        F.f = F.fexp.*F.f;
        if ~isempty(fired)      
        % update of synaptic filter with spikes
            F.f(fired) = F.f(fired) + F.fSpike(fired);
        end
    elseif strcmp(dynamicsModel,'dexp')
        F.f_r_AMPA = F.fexp_r_AMPA.*F.f_r_AMPA;
        F.f_d_AMPA = F.fexp_d_AMPA.*F.f_d_AMPA;
        F.f_r_NMDA = F.fexp_r_NMDA.*F.f_r_AMPA; 
        F.f_d_NMDA = F.fexp_d_NMDA.*F.f_d_NMDA;
        if ~isempty(fired)
            F.f_r_AMPA(fired) = F.f_r_AMPA(fired) + F.fSpike_r(fired);
            F.f_d_AMPA(fired) = F.f_d_AMPA(fired) + F.fSpike_d(fired);
            F.f_r_NMDA(fired) = F.f_r_NMDA(fired) + F.fSpike_r(fired);
            F.f_d_NMDA(fired) = F.f_d_NMDA(fired) + F.fSpike_d(fired);
        end
        F.f_AMPA = F.fScale_AMPA.*(F.f_r_AMPA + F.f_d_AMPA);
        F.f_NMDA = F.fScale_NMDA.*(F.f_r_NMDA + F.f_d_NMDA);
        F.f = 0.80*F.f_AMPA + 0.20*F.f_NMDA;
    end
end

end