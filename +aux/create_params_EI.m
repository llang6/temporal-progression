% Creates parameter files from choice of parameter set specified in Opt
% and saves it in .../DATA/Params.mat

function create_params_EI(paramsfile)
%----------
% NETWORKS
%----------
% Start and end of trial (in units of seconds)
Sim.t_Start = -1;
Sim.t_End = 3;
Sim.dt_step = 0.0001; % integration step (s)

%------------------
% NETWORK OPTIONS
%------------------
Network.clusters = {'EE','EI','IE','II','decision'}; % all weights are clustered
Network.clust = 'hom'; % homogeneous EE cluster size
% Network.clust = 'het'; % heterogeneous EE cluster size
Network.clust_std = 0.01; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
Network.clustEI = 'hom'; % EI clusters: homogeneous='hom', heterogeneous='het'
Network.clustIE = 'hom'; % IE clusters: homogeneous='hom', heterogeneous='het'
Network.clustII = 'hom'; % II clusters: homogeneous='hom', heterogeneous='het'
% Network.clust_syn = '';
N = 2000; % 2000 % network size
N_e = N*4/5; % (4/5)*N % exc neurons
N_i = N/5; % (1/5)*N % inh neurons
Scale = (1000/N)^(1/2);

% % global spontaneous firing rates (neeed to fix thresholds)
% ni_e = 2; % 3.6;   %6.6   % AB97: 3.0
% ni_i = 5; % 5.2;   %8.2   % AB97: 4.2
%------------------
% TIME CONSTANTS
%------------------
tau_arp = 0.005; % refractory period
tau_e = 0.02; % exc membrane time
tau_i = 0.02; % inh membrane time
tausyn_e = 0.005; % exc synaptic time
tausyn_i = 0.005; % inh synaptic time
tau_E_r_AMPA = 0.001;
tau_E_d_AMPA = 0.005;
tau_E_r_NMDA = 0.005;
tau_E_d_NMDA = 0.050;
tau_I_r_AMPA = 0.001;
tau_I_d_AMPA = 0.005;
tau_I_r_NMDA = 0.005;
tau_I_d_NMDA = 0.050;
%------------------
% SYNAPTIC WEIGHTS
%------------------
% syn weights are drawn from a gaussian distribution with std delta and
% mean Jab
Jee = Scale*0.020; % (E-->E) % Scale*0.02
Jii = Scale*0.120; % (I--|I) % Scale*0.12
Jie = Scale*0.020; % (E-->I) % Scale*0.02
Jei = Scale*0.060; % (I--|E) % Scale*0.06
%--------------------
% CLUSTER PARAMETERS
%--------------------
% delta = 0.01;
delta = 0; % SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Jplus = 13; % 13 % (E-->E) EE intra-cluster potentiation factor
Network.factorEI = 10; % 10 % (I--|E) EI intra-cluster potentiation factor
Network.factorIE = 8; % 8 % (E-->I) IE intra-cluster potentiation factor
Network.factorII = 5; % 5 % (I--|I) II intra-cluster potentiation factor

bgr = 0.1; % 0.1 % fraction of background excitatory neurons (unclustered)
Network.bgrI = 0.1; % 0.1 % fraction of background neurons
Ncluster = 100; % 100 % average # neurons per E cluster
p = round(N_e*(1-bgr)/Ncluster); % # of clusters
f = (1-bgr)/p;
Network.fI = (1-Network.bgrI)/p;       % 0.09
% gam = 0.5; % % parameter related to inter-cluster depression (see function aux.SynWeights)
% rho = 2.75;
if Jplus < 1 || Jplus > (2-(p+1)*(bgr+1))/(bgr-1)+1, fprintf('\nWARNING: Jplus is outside of the intended range. Unexpected behavior may result. Continue?\n\n'); pause; end
if Network.factorII < 1 || Network.factorII > (2-(p+1)*(bgr+1))/(Network.bgrI-1)+1, fprintf('\nWARNING: Network.factorII is outside of the intended range. Unexpected behavior may result. Continue?\n\n'); pause; end

%------------------
% THRESHOLDS
%------------------
theta_e = 1.42824; 
theta_i = 0.74342;
% theta_e = 1; % exc threshold potential
% theta_i = 1; % inh threshold potentials
% reset potentials
He = 0; %
Hi = 0; %
%--------------------------
% CONNECTIVITY PARAMETERS
%--------------------------
Cee = N_e*0.2; % # presynaptic neurons
Cie = N_e*0.5; %
Cii = N_i*0.5; %
Cei = N_i*0.5; %
%------------------
% EXTERNAL BIAS
%------------------
% external input parameters, eg: external current given by mu_e_ext = Cext*Jee_ext*ni_ext
Cext = (N_e)*0.2; % # presynaptic external neurons
Jie_ext = 0.8*Scale*0.0915; % external input synaptic strengths
Jee_ext = 0.8*Scale*0.1027; %
% EXTERNAL CURRENT
% default external currents
ni_ext = 5; % 7;
mu_E0 = Cext*Jee_ext*ni_ext;
mu_I0 = Cext*Jie_ext*ni_ext;
% random ext current, same in each trial of simulation
Mu = [mu_E0*(ones(N_e,1)+(0.1/2)*(2*rand([N_e,1])-1)); ...
      mu_I0*(ones(N_i,1)+(0.05/2)*(2*rand([N_i,1])-1))]; % bias
% Mu = [mu_E0*(ones(N_e,1)-(0.1/2)); ...
%       mu_I0*(ones(N_i,1)-(0.05/2))]; % bias
  
%----------------
% DEFAULT STIMULI
%----------------
% STIMULUS
Stimulus.input = 'Const'; % constant external current
scnt = 0;
% TASTE (specific stimulus)
scnt = scnt + 1;
feat(scnt).name = 'US'; % unconditioned stimulus (taste)
feat(scnt).interval = [0, Sim.t_End]; % stimulus interval
gain = 0.2; % stimulus value at 1 s
feat(scnt).gain = gain;
%feat(scnt).profile = @(t)t; % time course of stimulus, eg a linear ramp
feat(scnt).profile = @(t)1;
feat(scnt).selectivity = 'mixed';
feat(scnt).selective = rand(1,p) < 0.5; % US selective clusters
feat(scnt).connectivity = 0.5; % fraction of selective neurons within a selective cluster
% SUCROSE
scnt = scnt + 1;
feat(scnt).name = 'Sucrose'; 
feat(scnt).interval = [0 0.5]; % stimulus interval
gain_inc = 0.15; % %0.3 stimulus value at end for neurons that are excited in response
gain_dec = -0.11; % -0.1
feat(scnt).gain_inc = gain_inc;
feat(scnt).gain_dec = gain_dec;
feat(scnt).profile_inc = @(t)1;
feat(scnt).profile_dec = @(t)1;
feat(scnt).selectivity = 'mixed'; % not all excitatory clusters respond
% feat(scnt).selectiveFrac = 53/214; % 105/214 % fraction of clusters that respond
% feat(scnt).selective = rand(1,p) < 0.5; % T/F for selective clusters
% feat(scnt).fracInc = 25/53; % 35/105 % of those clusters that respond, what fraction are excited
feat(scnt).connectivity = 0.5; % fraction of selective neurons within a selective cluster
% QUININE
scnt = scnt + 1;
feat(scnt).name = 'Quinine'; 
feat(scnt).interval = [0 0.5]; % stimulus interval
gain_inc = 0.15; % stimulus value at end for neurons that are excited in response
gain_dec = -0.11;
feat(scnt).gain_inc = gain_inc;
feat(scnt).gain_dec = gain_dec;
feat(scnt).profile_inc = @(t)1;
feat(scnt).profile_dec = @(t)1;
feat(scnt).selectivity = 'mixed'; % not all excitatory clusters respond
% feat(scnt).selectiveFrac = 45/214; % 83/214 % fraction of clusters that respond
% feat(scnt).selective = rand(1,p) < 0.5; % T/F for selective clusters
% feat(scnt).fracInc = 21/45; % 28/83 % of those clusters that respond, what fraction are excited
feat(scnt).connectivity = 0.5; % fraction of selective neurons within a selective cluster
% MALTOSE
scnt = scnt + 1;
feat(scnt).name = 'Maltose'; 
feat(scnt).interval = [0 0.5]; % stimulus interval
gain_inc = 0.15; % stimulus value at end for neurons that are excited in response
gain_dec = -0.11;
feat(scnt).gain_inc = gain_inc;
feat(scnt).gain_dec = gain_dec;
feat(scnt).profile_inc = @(t)1;
feat(scnt).profile_dec = @(t)1;
feat(scnt).selectivity = 'mixed'; % not all excitatory clusters respond
% feat(scnt).selectiveFrac = 61/214; % 109/214 % fraction of clusters that respond
% feat(scnt).selective = rand(1,p) < 0.5; % T/F for selective clusters
% feat(scnt).fracInc = 30/61; % 39/109 % of those clusters that respond, what fraction are excited
feat(scnt).connectivity = 0.5; % fraction of selective neurons within a selective cluster
% OCTAACETATE
scnt = scnt + 1;
feat(scnt).name = 'Octaacetate'; 
feat(scnt).interval = [0 0.5]; % stimulus interval
gain_inc = 0.20; % stimulus value at end for neurons that are excited in response
gain_dec = -0.11;
feat(scnt).gain_inc = gain_inc;
feat(scnt).gain_dec = gain_dec;
feat(scnt).profile_inc = @(t)1;
feat(scnt).profile_dec = @(t)1;
feat(scnt).selectivity = 'mixed'; % not all excitatory clusters respond
% feat(scnt).selectiveFrac = 48/214; % 96/214 % fraction of clusters that respond
% feat(scnt).selective = rand(1,p) < 0.5; % T/F for selective clusters
% feat(scnt).fracInc = 22/48; % 32/196 % of those clusters that respond, what fraction are excited
feat(scnt).connectivity = 0.5; % fraction of selective neurons within a selective cluster

%------------------------------------------------------------------
% SELECTIVITY TO TASTES
% get from structure 'combo2.mat'
if ~exist('combo2','var'), load('combo2.mat'); end
%[selectivity, ~] = aux.getSelectiveClusters(combo2,p,2);
%selectivity = aux.getSelectiveClustersV2(combo2,p);
%[selectivity,ind_inc,ind_dec] = aux.getSelectiveClustersV3;
% [selectivity,ind_inc,ind_dec] = aux.getSelectiveClustersV4(combo2,true,0);
% feat(2).selective = selectivity(:,1);
% feat(2).ind_inc = ind_inc{1};
% feat(2).ind_dec = ind_dec{1};
% feat(3).selective = selectivity(:,2);
% feat(3).ind_inc = ind_inc{2};
% feat(3).ind_dec = ind_dec{2};
% feat(4).selective = selectivity(:,3);
% feat(4).ind_inc = ind_inc{3};
% feat(4).ind_dec = ind_dec{3};
% feat(5).selective = selectivity(:,4);
% feat(5).ind_inc = ind_inc{4};
% feat(5).ind_dec = ind_dec{4};
selectivity = aux.getSelectiveClustersV6(p);
feat(2).selective = selectivity(:,1);
feat(3).selective = selectivity(:,2);
feat(4).selective = selectivity(:,3);
feat(5).selective = selectivity(:,4);
%------------------------------------------------------------------
                
% ANTICIPATORY CUE
scnt = scnt+1;
feat(scnt).name = 'CSgauss'; % conditioned stimulus (cue) with "quenched" noise
feat(scnt).interval = [-0.5 Sim.t_End]; % stimulus interval
gain = 0.2; % SD of quenched noise across neurons
feat(scnt).gain = gain;
tau_cue = [0.5, 1]; % rise and decay time of double exp cue time course
feat(scnt).profile = @(t)(1/(tau_cue(2)-tau_cue(1)))*(exp(-t/tau_cue(2))-exp(-t/tau_cue(1))); % double exp profile time course
feat(scnt).selectivity = 'exc'; % targets exc neurons only
feat(scnt).selective = ones(1,p); % CS targets all clusters
feat(scnt).connectivity = 0.50; % fraction of neurons targeted within each cluster

Stimulus.feat = feat;

%%
%--------------------
% DECISION CLUSTERS
%--------------------
% Decision clusters send input only to each other
% gam = 1/(2 - f*(p + 1));
% Jminus = 1 - gam*f*(Jplus - 1);
% jDefault = Jminus*Jee*0.1; % default excitatory cluster input
% jHigh = Jplus*Jee*0.085; % preferential excitatory cluster input
% jdd_out = -Network.factorII*Jii; % inter cluster inhibition
% jdd_in = Jplus*Jee*0.45; % intra cluster potentiation
% pDefault = 0.2; % probability of connections
% pHigh = 0.5;
N_decision = 100; % number of neurons in each decision cluster
tau_decision = 0.02; % membrane time constant (for ref, tau_e and tau_i are 0.02)
tausyn_decision = 0.05; % synaptic time constant (for ref, tausyn_e and tausyn_i are 0.005)
%------------------------------------------
% Connection strengths and probabilities
%------------------------------------------
% Background (unclustered) E --> L, R
jBackground_e = 0;            pBackgroundE = 0.2;
% Background (unclustered) I --> L, R
jBackground_i = 0;            pBackgroundI = 0.5;
% Default excitatory cluster --> L, R
jDefault_e = 0.0036;          pDefaultE = 0.2;
% Default inhibitory cluster --> L, R
jDefault_i = -0.0258;         pDefaultI = 0.5;
% Selective clusters excited by stimulus --> correct decision
%jHigh_e = 0.1 * 1.5;          pPrefE = 0.2;
jHigh_e = 0.1 * 1.5 * 1.2;          pPrefE = 0.2;
% Selective clusters inhibited by stimulus --> correct decision
jHigh_i = -0.2582;            pPrefI = 0.5;
% L --> L and R --> R
%jdd_in = 0.1838 * 1.5 * 1.2;  pSelf = 0.2;
jdd_in = 0.1838 * 1.5 * 1.1;  pSelf = 0.2;
% L --| R and R --| L
%jdd_out = -0.2582 * 0.9;      pCross = 0.5; 
jdd_out = -0.2582 * 0.8;      pCross = 0.5;

%%
%------------------------------------------------------------------------
% PLOT PARAMETERS: ------------------------------------------------
%------------------------------------------------------------------------
% PLOTS
Sim.Plotf = 0;
Sim.plot_length = Sim.t_End - Sim.t_Start; % length of plot intervals
% indices of ensemble units to store
exc = randperm(N_e);
inh = N_e+1:N_e+N_i;
inh = inh(randperm(numel(inh)));
Sim.ind_p = [exc(1) inh(1)]; % choosing neuron index for membrane potential plot (one E and one I)
Sim.weights_save = 'off'; % save weight matrix: 'Yes'
extra = '';

save(paramsfile,'ni_ext','tau_arp','tau_i','tau_e','tau_decision','tau_E_r_AMPA',...
    'tau_E_d_AMPA','tau_E_r_NMDA','tau_E_d_NMDA','tau_I_r_AMPA','tau_I_d_AMPA','tau_I_r_NMDA',...
    'tau_I_d_NMDA','theta_e','theta_i','delta','f','Scale','Jee','Jii','Jie','Jei','Jee_ext','Jie_ext',...
    'Jplus','jBackground_e','jBackground_i','jDefault_e','jDefault_i','jHigh_e','jHigh_i',...
    'jdd_in','jdd_out','pBackgroundE','pBackgroundI','pDefaultE','pDefaultI','pPrefE','pPrefI',...
    'pSelf','pCross','He','Hi','N_e','N_i','N_decision','Cee','Cie','Cei','Cii','Cext','p',...
    'Sim','Network','Stimulus','tausyn_e','tausyn_i','tausyn_decision','extra','Mu','paramsfile');

fprintf('Network parameters saved in %s\n',paramsfile);