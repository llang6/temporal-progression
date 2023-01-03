function [SPIKE_DATA,INPUT_DATA] = simulation(NET,NEU,STIM,INDS,STIMULUS,GATE_SWITCH,GATE_CURVE,...
    S,IS_MODIFY_WEIGHTS,WARMUP,BIN_W,TIME_V,SILENCING,SAVE_INPUTS,V0)
%[SPIKE_DATA,INPUT_DATA] = simulation(NET,NEU,STIM,INDS,STIMULUS,GATE_SWITCH,GATE_CURVE,...
%    S,IS_MODIFY_WEIGHTS,WARMUP,BIN_W,TIME_V,SILENCING,SAVE_INPUTS,V0) 
%
%Main function for LIF network simulations (see lifNet for more information
%    about each component of the model)
%
%Simulates a LIF network over time, with (optional) external stimulus,
%    action gate, and optogenetic silencing inputs, then returns all
%    spiking data
%
% Input NET: structure, contains parameters that define the network
%     architecture
%
% Input NEU: structure, contains parameters that define single neuron
%     properties
%
% Input STIM: structure, contains all information about the stimulus except
%     which neurons receive it in this simulation
%
% Input INDS: structure, contains vectors of indices that tell the network
%     which neurons belong to which functional cluster
%
% Input STIMULUS: string, indicates which stimulus ('Sucrose', 'Maltose',
%     'Quinine', 'Octaacetate', or '' for none) is applied in this
%     simulation
%
% Input GATE_SWITCH: logical, indicates whether or not the action gate will
%     be applied in this simulation
%
% Input GATE_CURVE: vector, the time course of the action gate to apply in
%     this simulation
%
% Input S: matrix, the synaptic weight matrix
%
% Input IS_MODIFY_WEIGHTS: logical, indicates whether or not the time
%     constants of cue and action clusters will be modified in this
%     simulation (this was always done along with modifying the connection
%     weights, although weight modifications are done to S before passing
%     it in)
%
% Input WARMUP: integer, time (in ms) before simulation starts saving data
%
% Input BIN_W: integer, time (in ms) before simulation offloads data (for
%     efficiency, not a parameter of the network itself)
%
% Input TIME_V: vector, time (is ms) over which simulation runs (including
%     WARMUP)
%
% Input SILENCING: vector, time course of the simulated optogenetic
%     silencing to apply to inhibitory neurons in this simulation (pass []
%     here to turn silencing off)
%
% Input SAVE_INPUTS: logical, indicates whether cluster-averaged input
%     current will be saved over time in addition to the spiking data (this
%     will significantly slow things down)
%
% Input V0: vector, initial membrane potentials for neurons in this
%     simulation (pass [] here to use the default, which is
%     4*randn(NET.N,1))
%
% Output SPIKE_DATA: matrix, column 1 is all spike times (in ms) from all
%     neurons in the simulation, column 2 is the corresponding index of the
%     neuron that spiked at that time
%
% Output INPUT_DATA: matrix, each row is an E cluster and each column is a
%     timepoint, gives the cluster-averaged total input current over time
%
% -LL
%

% initialization
q = 1; fired = []; binfired = [];
if isempty(V0)
    V = 4*randn(NET.N,1); % initial membrane potential
else
    V = V0;
end
tauv = NEU.tauE + 0*V; % membrane time constants
tauv((1+NET.Ne):NET.N) = NEU.tauI;
refr = 0*V; 
IsynE = 0*V; % synaptic currents
IsynI = 0*V;
ActionGate = ones(length(V),1);
silencing = ones(length(V),1);
if SAVE_INPUTS, INPUT_DATA = zeros(NET.Q+1,length(TIME_V)); else, INPUT_DATA = []; end
% constant inputs
Iext = NET.Ieext + 0*V;
Iext((1+NET.Ne):NET.N) = NET.Iiext;
Iext0 = Iext;
% synaptic time constants
TAUSE = NEU.tausE + 0*V;
TAUSI = NEU.tausI + 0*V;
if IS_MODIFY_WEIGHTS
    % excitatory synaptic input to excitatory cue clusters
    TAUSE(INDS.cueLE) = NEU.tauCueE;
    TAUSE(INDS.cueRE) = NEU.tauCueE;
    % inhibitory synaptic input to excitatory cue clusters
    TAUSI(INDS.cueLE) = NEU.tauCueI;
    TAUSI(INDS.cueRE) = NEU.tauCueI;
    % excitatory synaptic input to excitatory action clusters
    TAUSE(INDS.actLE) = NEU.tauActionE;
    TAUSE(INDS.actRE) = NEU.tauActionE;
    % inhibitory synaptic input to excitatory action clusters
    TAUSI(INDS.actLE) = NEU.tauActionI;
    TAUSI(INDS.actRE) = NEU.tauActionI;
end
% evolve network over time
DT = diff(TIME_V(1:2));
for i = 1:length(TIME_V)
    % reset external input
    Iext = Iext0;
    % data chunking: every BIN_W, offload full bin of spike data and start over
    t = TIME_V(i);
    if t>0 && ~mod(t,BIN_W)
        SPIKE_DATA{1,q} = binfired;
        q = q + 1; 
        binfired = [];
    end
    % reset neurons that just spiked
    V(V>=NEU.Vspk) = NEU.Vr;
    % apply stimulus to external current
    if ~isempty(STIMULUS)
        if (t-WARMUP)/1000 >= STIM.tStart
            if strcmp(STIMULUS,'Sucrose'), IND = STIM.sucInd;
            elseif strcmp(STIMULUS,'Quinine'), IND = STIM.quiInd;
            elseif strcmp(STIMULUS,'Maltose'), IND = STIM.malInd;
            elseif strcmp(STIMULUS,'Octaacetate'), IND = STIM.octInd;
            end
            Iext(IND) = Iext(IND) + STIM.Curve(i)*Iext(IND);
        end
    end
    % update synaptic current
    IsynE = IsynE+(-IsynE*DT+sum(S(:,fired(fired<=NET.Ne | fired>NET.N)),2))./TAUSE;
    IsynI = IsynI+(-IsynI*DT+sum(S(:,fired(fired>NET.Ne & fired<=NET.N)),2))./TAUSI;
    % update membrane potential
    if ~isempty(SILENCING), silencing((1:NET.N)>NET.Ne) = SILENCING(i); end
    Itotal = Iext.*silencing + IsynE + IsynI;
    if SAVE_INPUTS
        for c = 1:NET.Q+1
            INPUT_DATA(c,i) = mean(Itotal(NET.indE(c):(NET.indE(c+1)-1)));
        end
    end
    if strcmp(GATE_SWITCH,'on') 
        ActionGate(INDS.actLE) = GATE_CURVE(i);
        ActionGate(INDS.actRE) = GATE_CURVE(i);
        V = V - V*DT./tauv + Itotal.*ActionGate*DT;
    else
        V = V - V*DT./tauv + Itotal*DT;
    end
    % absolute refractory period
    V(refr>0) = NEU.Vr;
    refr = max(-1,refr-DT);
    % boundary conditions (spiking)
    fired = find(V>=NEU.Vth);
    if ~isempty(fired) 
        V(fired) = NEU.Vspk; refr(fired) = NEU.taur; 
        binfired = [binfired;t-WARMUP+0*fired,fired];
    end
end
end
