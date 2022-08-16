%% miscellaneous
%
% This script contains additional code for generating other plots and
% carrying out small analyses included in the paper.
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Experiment/SessionInfo.mat
% /RawData/Simulation/simulationDataX.mat for X in 245:254
% /HMMData/Experiment/classifiedStates_exp.mat
% /HMMData/Experiment/classifiedStates_exp_shuff_circ.mat
% /HMMData/Experiment/classifiedStates_exp_shuff_swap.mat
% /HMMData/Experiment/HMM_expX.mat for X in 1:21 
% /HMMData/Simulation/classifiedStates_sim.mat
% /HMMData/Simulation/classifiedStates_sim_shuff_circ.mat
% /HMMData/Simulation/classifiedStates_sim_shuff_swap.mat
% /HMMData/Simulation/HMM_simX.mat for X in 245:254
% requires access to functions:
% distinguishable_colors (in +fun)
% getStateCell (in +fun)
% scoreTrials (in +fun)

%% count numbers of HMM states (for Supplementary Tables 1 and 2)

% ------------------------------------------------------------
fileName = 'classifiedStates_exp.mat';
% ------------------------------------------------------------

% setup and load data
if contains(fileName,'exp')
    fileName = [sprintf('%s/HMMData/Experiment/',pwd),fileName];
    sessions = setdiff(1:21,[5,7,17,19,20]);
elseif contains(fileName,'sim')
    fileName = [sprintf('%s/HMMData/Simulation/',pwd),fileName];
    sessions = 245:254;
end
stateData = load(fileName);
temp = fieldnames(stateData);
temp = temp{1};
stateData = stateData.(temp);

% Hidden and Decoded states
numHiddenStates = [];
numDecodedStates = [];
for i = sessions
    tf = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf,3);
    numHiddenStates = [numHiddenStates, numel(block)];
    numDecodedStates = [numDecodedStates, sum(cellfun(@(x)~strcmp(x,'Not decoded'),block))];
end
fprintf('\nHidden states');
fprintf('\n    Total: %i',sum(numHiddenStates));
fprintf('\n    Mean: %.1f',mean(numHiddenStates));
fprintf('\n    Median: %.1f',median(numHiddenStates));
fprintf('\n    Range: %i - %i\n',min(numHiddenStates),max(numHiddenStates));
fprintf('\nDecoded states');
fprintf('\n    Total: %i',sum(numDecodedStates));
fprintf('\n    Mean: %.1f',mean(numDecodedStates));
fprintf('\n    Median: %.1f',median(numDecodedStates));
fprintf('\n    Range: %i - %i\n',min(numDecodedStates),max(numDecodedStates));

% Decision-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Exclusive Decision-coding','Cue-coding','Action-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nDecision-coding states: %i over %i sessions\n',numStates,numSessions);

% Cue-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Cue-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nCue-coding states: %i over %i sessions\n',numStates,numSessions);

% Action-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Action-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nAction-coding states: %i over %i sessions\n',numStates,numSessions);

% Quality-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Exclusive Quality-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nQuality-coding states: %i over %i sessions\n',numStates,numSessions);

% Taste ID-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Taste ID-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nTaste ID-coding states: %i over %i sessions\n',numStates,numSessions);

% Dual-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Dual-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nDual-coding states: %i over %i sessions\n',numStates,numSessions);

% Non-coding states
numStates = 0;
sessions_found = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    tf2 = ismember(block,{'Non-coding'});
    numStates = numStates + sum(tf2);
    if sum(tf2)>0, sessions_found = [sessions_found, i]; end
end
numSessions = length(sessions_found);
fprintf('\nNon-coding states: %i over %i sessions\n',numStates,numSessions);

% How many sessions contain Quality-coding AND Decision-coding states?
numSessions = 0;
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Exclusive Decision-coding','Cue-coding','Action-coding'}))
        numSessions = numSessions + 1;
    end
end
fprintf('\nNumber of sessions with Quality-coding and Decision-coding states: %i\n',numSessions);

% How many sessions contain Cue-coding AND Action-coding states?
numSessions = 0;
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Cue-coding'})) && ...
            any(ismember(block,{'Action-coding'}))
        numSessions = numSessions + 1;
    end
end
fprintf('\nNumber of sessions with Cue-coding and Action-coding states: %i\n',numSessions);

% How many sessions contain Quality-, Cue-, AND Action-coding states?
numSessions = 0;
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Cue-coding'})) && ...
            any(ismember(block,{'Action-coding'}))
        numSessions = numSessions + 1;
    end
end
fprintf('\nNumber of sessions with Quality-, Cue-, and Action-coding states: %i\n',numSessions);

%% plot any experimental trial with HMM state overlay (for Figure 1b-d, left)

% Figure 1b (left) is Session 13, Trial 218 (state 8)
% Figure 1c (left) is Session 16, Trial 166 (state 2)
% Figure 1d (left) is Session 2, Trial 115 (state 2)

% ---------------------
session = 2;
trial = 115;
% ---------------------

classifiedData = load(sprintf('%s/HMMData/Experiment/classifiedStates_exp.mat',pwd));
classifiedData = classifiedData.classifiedStates_exp;
hmmData = load(sprintf('%s/HMMData/Experiment/HMM_exp%i.mat',pwd,session));
spikes = hmmData.spikes;
win_train = hmmData.win_train;
[~,nneurons] = size(spikes);
tDecision = diff(win_train(trial,:))-0.2;
spks = [];
for i = 1:nneurons
    times = spikes(trial,i).spk;
    spks = [spks ; times , i*ones(length(times),1)];
end
spks(:,1) = spks(:,1) + tDecision;

figure(1); clf; hold all;
colors = fun.distinguishable_colors(30);
plot([spks(:,1)';spks(:,1)'],[spks(:,2)'-0.2;spks(:,2)'+0.2],'k','linewidth',1.5);
seq = hmmData.res.hmm_postfit(trial).sequence;
for i = 1:size(seq,2)
    patch([seq(1,i),seq(2,i),seq(2,i),seq(1,i)]+tDecision,...
        [0.6,0.6,nneurons+0.4,nneurons+0.4],...
        colors(seq(4,i),:),'edgecolor','none','facealpha',0.25);
    tf = vertcat(classifiedData{2:end,1})==session & vertcat(classifiedData{2:end,2})==seq(4,i);
    if contains(classifiedData{[false;tf],3},'Quality')
        text(seq(1,i)+tDecision,nneurons+0.6,'Quality-coding','color','r');
    elseif contains(classifiedData{[false;tf],3},'Cue')
        text(seq(1,i)+tDecision,nneurons+0.6,'Cue-coding','color','c');
    elseif contains(classifiedData{[false;tf],3},'Action')
        text(seq(1,i)+tDecision,nneurons+0.6,'Action-coding','color','b');
    end 
end
prob = hmmData.res.hmm_results(trial).pStates;
t = linspace(-0.1,tDecision+0.1,size(prob,2));
for i = 1:size(prob,1)
    plot(t,prob(i,:)*(nneurons-0.2)+0.6,'color',colors(i,:),'linewidth',1.5);
end
xline(0,'k'); xline(tDecision,'k');
xlim([-0.1,tDecision+0.1]); ylim([0.6,nneurons+0.4]);
yticks(1:nneurons);
set(gca,'tickdir','out','color','none','box','off');
set(gcf,'renderer','painters');
%print(sprintf('Raster_S%iT%i',session,trial),'-dpdf');

%% plot any simulation trial (for Figure 3c and 4a)

% 3c (top) is Session 251, Baseline
% 3c (bottom) is Session 251, Trial 10
% 4a is Session 254, Trial 40

% -------------------------------------------------
session = 251;
trial = 'B'; % 'B' for baseline or # from 1-100
% -------------------------------------------------

data = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,session));
offsetX = -0.9;
width = 342;
height = 293;
if ischar(trial)
    % simulate a baseline trial (no stimulus) and plot it
    rng(99997);
    simulateAndPlotB(data,offsetX,width,height);
    xlabel('Time [s]','fontsize',10)
    ylabel('Neurons','fontsize',10);
    xlim([0,3]);
    yticks([]);
    title(sprintf('Session: %i, Baseline',session));
    set(gca,'fontsize',10,'Color','none','tickdir','out','box','off'); 
    set(gcf,'PaperPositionMode','Auto'); 
    set(gcf,'Renderer','painters');
    %print(sprintf('raster_sim_S%iB',session),'-dpng','-r0');
    %print(sprintf('raster_sim_S%iB',session),'-depsc');
else
    % plot trial raster plot
    firings_all = data.firings_all;
    firings = firings_all{trial};
    stimuli = data.stimuli;
    net = data.net;
    totalT = data.totalT;
    
    f = figure(2); clf; hold all;
    f.Position(3:4) = [width,height];
    rasterPlotAbbrev(firings,data.inds,data.clustersWithRoles,totalT,net,7,true,'tick',offsetX);
    xlabel('Time [s]','fontsize',10)
    ylabel('Neurons','fontsize',10);
    xlim([0,3]);
    yticks([]);
    title(sprintf('Session: %i, Trial: %i (%s)',session,trial,stimuli{trial}));
    set(gca,'fontsize',10,'Color','none','tickdir','out','box','off'); 
    set(gcf,'PaperPositionMode','Auto'); 
    set(gcf,'Renderer','painters');
    %print(sprintf('raster_sim_S%iT%i',session,trial),'-dpng','-r0');
    %print(sprintf('raster_sim_S%iT%i',session,trial),'-depsc');
end

%% plot any simulation trial with HMM state overlay (for Figure 4b and c-e ,left)

% 4b is Session 254, Trial 40
% 4c (left) is Session 247, Trial 58 (state 20)
% 4d (left) is Session 250, Trial 61 (state 25)
% 4e (left) is Session 249, Trial 72 (state 6)

% -----------------
session = 254;
trial = 40;
% -----------------

classifiedData = load(sprintf('%s/HMMData/Simulation/classifiedStates_sim.mat',pwd));
classifiedData = classifiedData.classifiedStates_sim;
hmmData = load(sprintf('%s/HMMData/Simulation/HMM_sim%i.mat',pwd,session));
networkData = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,session));
spikes = hmmData.spikes;
[~,nneurons] = size(spikes);
cOrder = [networkData.clustersWithRoles,setdiff(1:14,networkData.clustersWithRoles)];
plotInd = NaN(1,14); for i = 1:nneurons, plotInd(i) = find(cOrder==i); end
spks = [];
for i = 1:nneurons
    times = spikes(trial,i).spk;
    spks = [spks ; times , plotInd(i)*ones(length(times),1)];
end

figure(3); clf; hold all;
colors = fun.distinguishable_colors(30);
plot([spks(:,1)';spks(:,1)'],[spks(:,2)'-0.2;spks(:,2)'+0.2],'k','linewidth',1.5);
seq = hmmData.res.hmm_postfit(trial).sequence;
for i = 1:size(seq,2)
    patch([seq(1,i),seq(2,i),seq(2,i),seq(1,i)],...
        [0.6,0.6,nneurons+0.4,nneurons+0.4],...
        colors(seq(4,i),:),'edgecolor','none','facealpha',0.25);
    tf = vertcat(classifiedData{2:end,1})==session & vertcat(classifiedData{2:end,2})==seq(4,i);
    if contains(classifiedData{[false;tf],3},'Quality')
        text(seq(1,i),nneurons+0.9,'Quality-coding','color','r');
    elseif contains(classifiedData{[false;tf],3},'Cue')
        text(seq(1,i),nneurons+0.9,'Cue-coding','color','c');
    elseif contains(classifiedData{[false;tf],3},'Action')
        text(seq(1,i),nneurons+0.9,'Action-coding','color','b');
    end 
end
prob = hmmData.res.hmm_results(trial).pStates;
t = linspace(0,3,size(prob,2));
for i = 1:size(prob,1)
    plot(t,prob(i,:)*(nneurons-0.2)+0.6,'color',colors(i,:),'linewidth',1.5);
end
xlim([0,3]); ylim([0.6,nneurons+0.4]);
yticks(1:nneurons);
set(gca,'tickdir','out','color','none','box','off');
set(gcf,'renderer','painters');
%print(sprintf('raster_sim_hmm_S%iT%i',session,trial),'-dpdf');

%% plot any state's occurrence frequencies given trial types for any experimental session (for Figure 1b-d, right)

% 1b (right) is Session 13, State 8
% 1c (right) is Session 16, State 2
% 1d (right) is Session 2, State 2

% ----------------
session = 13;
state = 8;
% ----------------

hmmData = load(sprintf('%s/HMMData/Experiment/HMM_exp%i.mat',pwd,session));
behaviorData = load(sprintf('%s/RawData/Experiment/SessionInfo.mat',pwd));
stimuli = {SessionInfo(session).TrialEvents.TrialType(:).TasteID};
scoredTrials = [SessionInfo(session).TrialEvents.TrialType(:).wasCorrect];
myCell = fun.getStateCell(hmmData,stimuli,scoredTrials,true);

f = figure(4); clf;
f.Position(3:4) = [332,136];
correct = [sum([myCell{[1,2],1}]==state)/sum([myCell{[1,2],2}]) sum([myCell{[3,4],1}]==state)/sum([myCell{[3,4],2}])];
incorrect = [sum([myCell{[1,2],3}]==state)/sum([myCell{[1,2],4}]) sum([myCell{[3,4],3}]==state)/sum([myCell{[3,4],4}])];
bar([correct' incorrect']);
set(gca,'XTickLabel',{'Cue Left','Cue Right'},'fontsize',10);
set(gca,'Color','none','TickLength',[0 0],'box','off');
legend('Correct Trials','Incorrect Trials','location','northeastoutside','fontsize',10);
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'PaperPositionMode','auto')

f = figure(5); clf;
f.Position(3:4) = [212,136];
correct = [sum([myCell{[1,3],1}]==state)/sum([myCell{[1,3],2}]) sum([myCell{[2,4],1}]==state)/sum([myCell{[2,4],2}])];
bar(correct');
set(gca,'XTickLabel',{'Sweet','Bitter'},'fontsize',10);
set(gca,'Color','none','TickLength',[0 0],'box','off');
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'PaperPositionMode','auto')

%% plot any state's occurrence frequencies given trial types for any simulation session (for Figure 4c-e, right)

% 4c (right) is Session 247, state 20
% 4d (right) is Session 250, state 25
% 4e (right) is Session 249, state 6

% ----------------
session = 249;
state = 6;
% ----------------

hmmData = load(sprintf('%s/HMMData/Simulation/HMM_sim%i.mat',pwd,session));
networkData = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,session));
classificationData = load(sprintf('%s/HMMData/Simulation/classifiedStates_sim.mat',pwd));
scoredTrials = fun.scoreTrials(networkData);
myCell = fun.getStateCell(hmmData,networkData.stimuli,scoredTrials,false);

f = figure(6); clf;
f.Position(3:4) = [332,136];
correct = [sum([myCell{[1,2],1}]==state)/sum([myCell{[1,2],2}]) sum([myCell{[3,4],1}]==state)/sum([myCell{[3,4],2}])];
incorrect = [sum([myCell{[1,2],3}]==state)/sum([myCell{[1,2],4}]) sum([myCell{[3,4],3}]==state)/sum([myCell{[3,4],4}])];
bar([correct' incorrect']);
set(gca,'XTickLabel',{'Cue Left','Cue Right'},'fontsize',10);
set(gca,'Color','none','TickLength',[0 0],'box','off');
legend('Correct Trials','Incorrect Trials','location','northeastoutside','fontsize',10);
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'PaperPositionMode','auto')

f = figure(7); clf;
f.Position(3:4) = [212,136];
correct = [sum([myCell{[1,3],1}]==state)/sum([myCell{[1,3],2}]) sum([myCell{[2,4],1}]==state)/sum([myCell{[2,4],2}])];
bar(correct');
set(gca,'XTickLabel',{'Sweet','Bitter'},'fontsize',10);
set(gca,'Color','none','TickLength',[0 0],'box','off');
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'PaperPositionMode','auto')

%% local function definitions
function simulateAndPlotB(DATA,X_OFFSET,WIDTH,HEIGHT)
    net = DATA.net;
    neu = DATA.neu;
    totalT = DATA.totalT;
    warmup = DATA.warmup;
    dt = DATA.dt;
    BinW = DATA.BinW;
    timev = 0:dt:totalT;
    S = DATA.S;
    gate = DATA.gate;
    inds = DATA.inds;
    gate.Curve = getActionGate(gate.StartMean,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    
    spkdata = simulation(net,neu,[],inds,'',gate,S,true,warmup,BinW,timev);
    
    f = figure(2); clf; hold all;
    f.Position(3:4) = [WIDTH,HEIGHT];
    firings = []; for i = 1:length(spkdata), firings = [firings; spkdata{i}]; end
    rasterPlotAbbrev(firings,inds,DATA.clustersWithRoles,totalT,net,7,true,'tick',X_OFFSET);
end

function CURVE = getActionGate(T_START,DURATION,FLOOR,CEILING,WARMUP,TOTAL_T,DT)
% create a ramp function that is FLOOR before time T_START [s], then ramps up
%     to CEILING over the course of DURATION [s]
% WARMUP, TOTAL_T, and DT are the time parameters we need to form the trial
%     time vector. Note that these are in [ms] and the ramp parameters are
%     converted to [ms] inside this code
timeV = 0:DT:(WARMUP+TOTAL_T);
CURVE = zeros(1,length(timeV));
CURVE(timeV>=WARMUP+1000*T_START) = min(1/(1000*DURATION)*(0:DT:(TOTAL_T-1000*T_START)),1);
CURVE = (CEILING-FLOOR)*CURVE + FLOOR;
end

function SPIKE_DATA = simulation(NET,NEU,STIM,INDS,STIMULUS,GATE,S,IS_MODIFY_WEIGHTS,WARMUP,BIN_W,TIME_V)
% initialization
q = 1; fired = []; binfired = [];
V = 0 + 4*randn(NET.N,1); % initial membrane potential
tauv = NEU.tauE + 0*V; % membrane time constants
tauv((1+NET.Ne):NET.N) = NEU.tauI;
refr = 0*V; 
IsynE = 0*V; % synaptic currents
IsynI = 0*V;
ActionGate = ones(length(V),1);
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
    % reset NEUrons that just spiked
    V(V>=NEU.Vspk) = NEU.Vr;
    % apply STIMULUS to external current
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
    if strcmp(GATE.Switch,'on') 
        ActionGate(INDS.actLE) = GATE.Curve(i);
        ActionGate(INDS.actRE) = GATE.Curve(i);
        V = V-V*DT./tauv+(Iext+IsynE+IsynI).*ActionGate*DT;
    else
        V = V-V*DT./tauv+(Iext+IsynE+IsynI)*DT;
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

function rasterPlotAbbrev(SPIKE_DATA,INDS,CLUSTERS_WITH_ROLES,TOTAL_T,NET,DOWNSAMPLE,APPROX_RATIO,DOT_OR_TICK,X_OFFSET)
downsampleEC = DOWNSAMPLE;
anchor = 0;
if APPROX_RATIO
    downsampleIC = max(round(downsampleEC/250*71),1);
    downsampleEB = 2*downsampleEC;
else
    downsampleIC = downsampleEC;
    downsampleEB = downsampleEC;
end 
% plot E taste clusters
ind_net = [INDS.sucE(randperm(length(INDS.sucE),downsampleEC)),...
    INDS.quiE(randperm(length(INDS.quiE),downsampleEC)),...
    INDS.malE(randperm(length(INDS.malE),downsampleEC)),...
    INDS.octE(randperm(length(INDS.octE),downsampleEC))];
ind_remap = (1:4*downsampleEC)'+anchor;
anchor = anchor + 4*downsampleEC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'k.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'k');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2,'Suc','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+downsampleEC,'Qui','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+2*downsampleEC,'Mal','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+3*downsampleEC,'Oct','VerticalAlignment','middle','fontsize',10);
% plot E cue left and action left clusters
ind_net = [INDS.cueLE(randperm(length(INDS.cueLE),downsampleEC)),...
    INDS.actLE(randperm(length(INDS.actLE),downsampleEC))];
ind_remap = [(1:downsampleEC)';((2*downsampleEC+1):3*downsampleEC)']+anchor;
anchor = anchor + downsampleEC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'g.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'g');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+4*downsampleEC,'Cue L','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+6*downsampleEC,'Act L','VerticalAlignment','middle','fontsize',10);
% plot E cue right and action right clusters
ind_net = [INDS.cueRE(randperm(length(INDS.cueRE),downsampleEC)),...
    INDS.actRE(randperm(length(INDS.actRE),downsampleEC))];
ind_remap = [(1:downsampleEC)';((2*downsampleEC+1):3*downsampleEC)']+anchor;
anchor = anchor + 3*downsampleEC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'b.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'b');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+5*downsampleEC,'Cue R','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+7*downsampleEC,'Act R','VerticalAlignment','middle','fontsize',10);
% plot remaining E clusters
for clust = setdiff(1:NET.Q,CLUSTERS_WITH_ROLES)
    ind_net = NET.indE(clust):(NET.indE(clust+1)-1); 
    ind_net = ind_net(randperm(length(ind_net),downsampleEC));
    ind_remap = (1:downsampleEC)'+anchor;
    anchor = anchor + downsampleEC;
    times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
    ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
    if strcmp(DOT_OR_TICK,'dot')
        plot(times,ind_plot,'k.','markersize',0.1);
    elseif strcmp(DOT_OR_TICK,'tick')
        times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
        plot([times;times],[ind_plot-0.4;ind_plot+0.4],'k');
    else
        error('Invalid option for DOT_OR_TICK. Use lowercase letters');
    end
end
% plot background E neurons
ind_net = NET.indE(NET.Q+1):(NET.indE(NET.Q+2)-1);
ind_net = ind_net(randperm(length(ind_net),downsampleEB));
ind_remap = (1:downsampleEB)'+anchor;
anchor = anchor + downsampleEB;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'k.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'k');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot I taste clusters
ind_net = [INDS.sucI(randperm(length(INDS.sucI),downsampleIC)),...
    INDS.quiI(randperm(length(INDS.quiI),downsampleIC)),...
    INDS.malI(randperm(length(INDS.malI),downsampleIC)),...
    INDS.octI(randperm(length(INDS.octI),downsampleIC))];
ind_remap = (1:4*downsampleIC)'+anchor;
anchor = anchor + 4*downsampleIC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'r.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot I cue left and action left clusters
ind_net = [INDS.cueLI(randperm(length(INDS.cueLI),downsampleIC)),...
    INDS.actLI(randperm(length(INDS.actLI),downsampleIC))];
ind_remap = [(1:downsampleIC)';((2*downsampleIC+1):3*downsampleIC)']+anchor;
anchor = anchor + downsampleIC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'r.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot I cue right and action right clusters
ind_net = [INDS.cueRI(randperm(length(INDS.cueRI),downsampleIC)),...
    INDS.actRI(randperm(length(INDS.actRI),downsampleIC))];
ind_remap = [(1:downsampleIC)';((2*downsampleIC+1):3*downsampleIC)']+anchor;
anchor = anchor + 3*downsampleIC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'r.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot remaining I clusters
for clust = setdiff(1:NET.Q,CLUSTERS_WITH_ROLES)
    ind_net = NET.indI(clust):(NET.indI(clust+1)-1); 
    ind_net = ind_net(randperm(length(ind_net),downsampleIC));
    ind_remap = (1:downsampleIC)'+anchor;
    anchor = anchor + downsampleIC;
    times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
    ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
    if strcmp(DOT_OR_TICK,'dot')
        plot(times,ind_plot,'r.','markersize',0.1);
    elseif strcmp(DOT_OR_TICK,'tick')
        times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
        plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
    else
        error('Invalid option for DOT_OR_TICK. Use lowercase letters');
    end
end
% visual formatting
yline(NET.Q*downsampleEC+downsampleEB+0.5,'r');
yline(NET.Q*(downsampleEC+downsampleIC)+downsampleEB+0.5,'r');
%set(gca,'fontsize',18,'tickdir','out','color','none','box','off');
%yticks([]);
xlim([0, TOTAL_T/1000]); ylim([0,NET.Q*(downsampleEC+downsampleIC)+downsampleEB+1]);
end
