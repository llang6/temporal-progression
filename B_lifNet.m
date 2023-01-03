%% lifNet
%
% Simulate a network of Leaky Integrate-and-Fire neurons that has been
% adapted to model performance of the perceptual decision-making task from
% Vincis et al. (2020).
%
% The basic options allow you to turn on or off the key model features as
% well as automatic plotting and saving.
%
% All other parameters can be altered in the 'parameters' section.
%
% Simulating 4 trials takes about 0.5-1 minute.
%
% Generated files were saved in '/RawData/Simulation'
%
% -LL
%

%% dependencies
% None
% This script relies on locally-defined functions only so that it can be
% self-contained. However, many of the local functions also exist in +fun,
% and calls to them could be prefixed by 'fun.' without changing the
% behavior

%% options
%close all; clearvars; clc; % start fresh
rng(9995); % make simulation results repeatable
isOverlapTastes = true;
isModifyWeights = true; % this also controls modified tau_syns for cue/action clusters
isActionGate = true;
isSilencing = false;
isPlot = true;
saveInputs = false; % saves total current input over time for clusters (will slow things down)
numTrials = 4;

%% parameters
% ----neuron----
neu.VL = 0; % resting potential, [mV]
neu.Vth = 20; % spiking threshold, [mV]
neu.Vr = 0; % reset voltage, [mV]
neu.C = 1; % membrane capacitance, [pF]
neu.tauE = 20; % E membrane time constant, [ms]
neu.tauI = 10; % I membrane time constant, [ms]
neu.tausE = 3; % time constant for E synaptic input, [ms]
neu.tausI = 2; % time constant for I synaptic input, [ms]
neu.taur = 5; % absolute refractory period, [ms]
neu.IthE = (neu.Vth-neu.VL)/neu.tauE*neu.C; % threshold current for E neurons, [pA]
neu.IthI = (neu.Vth-neu.VL)/neu.tauI*neu.C; % threshold current for I neurons, [pA]
neu.Vspk = 45; % overshoot of action potential, [mV] 

% ----network----
net.Ne = 4000; % number of E neurons
net.Ni = 994; % number of I neurons
net.f = 0.875;
net.N = net.Ne+net.Ni; % number of total neurons
net.J = 1/sqrt(net.N); % global scaling factor
net.Pei = 0.5; % mean J(inh --> exc) connectivity (probability of a connection)
net.Pee = 0.2; % mean J(exc --> exc) connectivity 
net.Pii = 0.5; % mean J(inh --> inh) connectivity
net.Pie = 0.5; % mean J(exc --> inh) connectivity
net.Jee = 0.16; % E->E synaptic weight, [pA]
net.Jei = -0.63; % I->E synaptic weight, [pA]
net.Jie = 0.23; % E->I synaptic weight, [pA]
net.Jii = -1.2; % I->I synaptic weight, [pA]
net.Ieext = 2.05; % external current for E neurons
net.Iiext = 2.16; % external current for I neurons
net.Q = 14; % number of clusters

% ----clusters and overlaps----
clusterRoles = {'S','Q','M','O','CL','CR','AL','AR'};
tasteOverlapFraction = 0.25; % fraction of neurons in each taste cluster selected for overlapping subcluster
tasteOverlapFactor = 0.75; % scale factor, between 0 and 1, to go from Jplus to Jplusminus
                           % higher factor ==> less contrast between overlapping and nonoverlapping neurons 

% ----connection modifications----
%     'Pot': potentiation (increased excitation) 
%     'Dep': depression (decreased excitation) 
%     'Inh': inhibition (increased inhibition)
% stimulus to cue
probPot_S2C = 0.6; % fraction of cue cluster neurons that receives boosted input from relevant stimulus cluster
factPot_S2C = 1.850; % cue clusters receive boosted input from relevant stimulus clusters
% cue to cue
probPot_C2C = 0.6;
probInh_C2C = 0.5;
factPot_C2C = 1.050; % modify self-excitation to promote persistant activity
factInh_C2C = 1.400; % 
% cue to action
probPot_C2Acor = 0.5; % fraction of action cluster neurons that receive boosted input from relevant cue cluster
probPot_C2Ainc = 0.5;
factPot_C2Acor = 2.750; % action clusters receive boosted input from relevant cue cluster...
factPot_C2Ainc = 2.600; % ...and also from irrelevant one, but less (more similar to above ==> more errors)
% action to cue
probInh_A2C = 0.5;
factInh_A2C = 1.400; % 
% action to action
probPot_A2A = 0.5;
probInh_A2A = 0.5;
factPot_A2A = 1.070; % modify self-excitation to promote persistant activity
factInh_A2A = 3.000; % 

% ----synaptic time constants----
neu.tauCueE = neu.tausE*2.300; % synaptic time constant for cue clusters' E input, [ms]
neu.tauCueI = neu.tausI*1.300; % synaptic time constant for cue clusters' I input, [ms]
neu.tauActionE = neu.tausE*3.300; % synaptic time constant for action clusters' E input, [ms]
neu.tauActionI = neu.tausI*3.250; % synaptic time constant for action clusters' I input, [ms]

% ----time----
warmup = 0; % simulate for this long before considering the trial to have begun, [ms]
totalT = 3000; % trial time, [ms]
dt = 0.05; % numerical integration time step (finest level of temporal resolution), [ms]
BinW = 10; % simulated spike data is collected in bins of this length (not important for simulation itself), [ms]

% ----stimulus----
stim.Profile = 'alpha'; % stimulus shape % 'square', 'alpha'
stim.tStart = 0; % stimulus begins this long after trial start, [s]
stim.Gain = 0.60; % maximum fractional increase in external current
stim.Duration = 0.5; % stimulus duration, [s] (relevant only if Profile is 'square')
stim.tauRise = 0.150; % rise time constant, [s] (relevant only if Profile is 'alpha')
stim.tauFall = 0.160; % fall time constant, [s] (relevant only if Profile is 'alpha')
stim.frac = 0.5; % fraction of responsive neurons in stimulus clusters

% ----action gate----
gate.Floor = 0; % minimum ramp value
gate.Ceiling = 1; % maximum ramp value
gate.StartMean = 1.5; % average start time of ramp, [s]
gate.Noise = 0.4; % half-width of uniform distribution from which start times are sampled, [s]
gate.Duration = 0.5; % duration of the ramp, [s]

% ----silencing----
silence.start = 1.25; % relative to stimulus onset, [s]
silence.duration = 0.5; % [s]
silence.strength = 100; % [% increase]

% ----performance scoring----
binSize = 0.05; % bin size for firing rate calculations, [s]
stateThresh = 0.4; % cluster is considered 'on' when its firing rate >= this fraction of its max firing rate...
stateMinTime = 0.05; % ...for at least this long, [s]
bound1 = 0.5; % time after stimulus when we start checking for action cluster activation, [s]
bound2 = 3.0; % time after stimulus when we stop checking for action cluster activation, [s]

%% setup 
timerStart_total = tic;
rngAtStart = rng;
homeDir = pwd; addpath(homeDir);

%% divide network into Q clusters (includes the background E neurons)
if mod(net.Ne*net.f,net.Q)
    error('Choose the proper number of clusters.')
end
net.indE = (1:floor(net.Ne*net.f/net.Q):net.Ne*net.f); 
if net.f == 1
    net.indE = [net.indE,net.Ne+1];
else
    net.indE = [net.indE,net.Ne*net.f+1,net.Ne+1];
end
net.Necl = diff(net.indE); % number of neurons in each E cluster
net.indI = (1:floor(net.Ni/net.Q):net.Ni); 
net.indI(end+1) = net.Ni+1; net.indI = net.indI+net.Ne;
net.Nicl = diff(net.indI); % number of neurons in each I cluster

%% synaptic matrix
fprintf('\nConstructing synaptic matrix ... ');
% construct the matrix
Jpe = 12;
Jme = max(0,(net.Q-net.f*Jpe)/(net.Q-net.f));
RJ = 0.4; % proportional factor between E and I clusters
Jpi = 1 + RJ*(Jpe-1); % same factor is also used for private I->E and E->I connections
Jmi = max(0,(net.Q-Jpi)/(net.Q-1));
S = [net.Jee*Jme*ones(net.Ne,net.Ne),net.Jei*Jmi*ones(net.Ne,net.Ni);...
     net.Jie*Jmi*ones(net.Ni,net.Ne),net.Jii*Jmi*ones(net.Ni,net.Ni)];
for i = 1:net.Q
    % within same E cluster
    S(net.indE(i):(net.indE(i+1)-1),net.indE(i):(net.indE(i+1)-1)) = net.Jee*Jpe;
    % within same I cluster
    S(net.indI(i):(net.indI(i+1)-1),net.indI(i):(net.indI(i+1)-1)) = net.Jii*Jpi;
    % E cluster -> its private I cluster
    S(net.indI(i):(net.indI(i+1)-1),net.indE(i):(net.indE(i+1)-1)) = net.Jie*Jpi;
    % I cluster -> its private E cluster
    S(net.indE(i):(net.indE(i+1)-1),net.indI(i):(net.indI(i+1)-1)) = net.Jei*Jpi;
end
% remove connections according to connectivity, taking overlaps into account 
clustersWithRoles = randperm(net.Q,length(clusterRoles));
overlappingClusters = clustersWithRoles(cellfun(@(x)ismember(x,{'S','Q','M','O'}),clusterRoles)); % overlap taste clusters
overlappingIndicesE = struct(); overlappingIndicesI = struct();
for i = 1:length(overlappingClusters)
    clust = overlappingClusters(i);
    overlappingIndicesE(i).cluster = clust;
    overlappingIndicesI(i).cluster = clust;
    overlappingIndicesE(i).all = net.indE(clust):net.indE(clust+1)-1;
    overlappingIndicesI(i).all = net.indI(clust):net.indI(clust+1)-1;
    poolE = overlappingIndicesE(i).all;
    poolI = overlappingIndicesI(i).all;
    overlappingIndicesE(i).overlap = poolE(randperm(length(poolE),round(tasteOverlapFraction*length(poolE))));
    overlappingIndicesI(i).overlap = poolI(randperm(length(poolI),round(tasteOverlapFraction*length(poolI))));
    overlappingIndicesE(i).nonoverlap = setdiff(poolE,overlappingIndicesE(i).overlap);
    overlappingIndicesI(i).nonoverlap = setdiff(poolI,overlappingIndicesI(i).overlap);
end
P = getP(net,overlappingIndicesE,overlappingIndicesI); % value of 0 in P marks corresponding connection in S for removal
S(~P) = 0;
fprintf('Synaptic matrix constructed.\n');

%% modify the synaptic matrix
% get the clusters with roles and their indices
clustSuc = clustersWithRoles(cellfun(@(x)strcmp(x,'S'),clusterRoles));
clustMal = clustersWithRoles(cellfun(@(x)strcmp(x,'M'),clusterRoles));
clustQui = clustersWithRoles(cellfun(@(x)strcmp(x,'Q'),clusterRoles));
clustOct = clustersWithRoles(cellfun(@(x)strcmp(x,'O'),clusterRoles));
clustCueL = clustersWithRoles(cellfun(@(x)strcmp(x,'CL'),clusterRoles));
clustCueR = clustersWithRoles(cellfun(@(x)strcmp(x,'CR'),clusterRoles));
clustActL = clustersWithRoles(cellfun(@(x)strcmp(x,'AL'),clusterRoles));
clustActR = clustersWithRoles(cellfun(@(x)strcmp(x,'AR'),clusterRoles));
inds.sucE = net.indE(clustSuc):(net.indE(clustSuc+1)-1); 
inds.sucI = net.indI(clustSuc):(net.indI(clustSuc+1)-1);
inds.malE = net.indE(clustMal):(net.indE(clustMal+1)-1); 
inds.malI = net.indI(clustMal):(net.indI(clustMal+1)-1);
inds.quiE = net.indE(clustQui):(net.indE(clustQui+1)-1); 
inds.quiI = net.indI(clustQui):(net.indI(clustQui+1)-1);
inds.octE = net.indE(clustOct):(net.indE(clustOct+1)-1); 
inds.octI = net.indI(clustOct):(net.indI(clustOct+1)-1);
inds.cueLE = net.indE(clustCueL):(net.indE(clustCueL+1)-1); 
inds.cueLI = net.indI(clustCueL):(net.indI(clustCueL+1)-1);
inds.cueRE = net.indE(clustCueR):(net.indE(clustCueR+1)-1); 
inds.cueRI = net.indI(clustCueR):(net.indI(clustCueR+1)-1);
inds.actLE = net.indE(clustActL):(net.indE(clustActL+1)-1); 
inds.actLI = net.indI(clustActL):(net.indI(clustActL+1)-1);
inds.actRE = net.indE(clustActR):(net.indE(clustActR+1)-1); 
inds.actRI = net.indI(clustActR):(net.indI(clustActR+1)-1);

% ----OVERLAPS----
if isOverlapTastes
    S = overlap(clustSuc,clustMal,...
                tasteOverlapFactor,S,overlappingIndicesE,overlappingIndicesI,net,Jme,Jpe,Jmi,Jpi);
    S = overlap(clustQui,clustOct,...
                tasteOverlapFactor,S,overlappingIndicesE,overlappingIndicesI,net,Jme,Jpe,Jmi,Jpi);
    fprintf('\nTaste clusters overlapped.\n');
end

% ----BIASED WEIGHTS----
if isModifyWeights
    % Stim L ==> Cue L, Stim R ==> Cue R (E->E pot Stim->Cue)
    S = modifyWeights(inds.cueLE,[inds.sucE,inds.quiE],factPot_S2C,probPot_S2C,S);
    S = modifyWeights(inds.cueRE,[inds.malE,inds.octE],factPot_S2C,probPot_S2C,S);
    % Cue L ==> Cue L, Cue R ==> Cue R (E->E pot Cue->Cue)
    S = modifyWeights(inds.cueLE,inds.cueLE,factPot_C2C,probPot_C2C,S);
    S = modifyWeights(inds.cueRE,inds.cueRE,factPot_C2C,probPot_C2C,S);
    % Cue L ==| Cue R, Cue R ==| Cue L (I-|E inh Cue->Cue)
    S = modifyWeights(inds.cueRE,inds.cueLI,factInh_C2C,probInh_C2C,S);
    S = modifyWeights(inds.cueLE,inds.cueRI,factInh_C2C,probInh_C2C,S);
    % Cue L ==> Action L, Cue R ==> Action R (E->E pot Cue->ActCor)
    S = modifyWeights(inds.actLE,inds.cueLE,factPot_C2Acor,probPot_C2Acor,S);
    S = modifyWeights(inds.actRE,inds.cueRE,factPot_C2Acor,probPot_C2Acor,S);
    % Cue L ==> Action R, Cue R ==> Action L (E->E pot Cue->ActInc)
    S = modifyWeights(inds.actRE,inds.cueLE,factPot_C2Ainc,probPot_C2Ainc,S);
    S = modifyWeights(inds.actLE,inds.cueRE,factPot_C2Ainc,probPot_C2Ainc,S);
    % Action L, Action R ==| Cue L, Cue R (I-|E inh Act->Cue)
    S = modifyWeights([inds.cueLE,inds.cueRE],[inds.actLI,inds.actRI],factInh_A2C,probInh_A2C,S);
    % Action L ==> Action L, Action R ==> Action R (E->E pot Act->Act)
    S = modifyWeights(inds.actLE,inds.actLE,factPot_A2A,probPot_A2A,S);
    S = modifyWeights(inds.actRE,inds.actRE,factPot_A2A,probPot_A2A,S);    
    % Action L ==| Action R, Action R ==| Action L (I-|E inh Act->Act)
    S = modifyWeights(inds.actRE,inds.actLI,factInh_A2A,probInh_A2A,S);
    S = modifyWeights(inds.actLE,inds.actRI,factInh_A2A,probInh_A2A,S);
    fprintf('\nSynaptic weights modified.\n');
end

%% re-index for plotting
% reMapCell is used only for the weight summary plot
% indexReMap can be used to get the 'real' index of the neuron in terms of
% its plotted index
% indexReMap_inv can be used to get the plotted index of the neuron in
% terms of its 'real' index
[reMapCell,indexReMap,indexReMap_inv] = reMapIndex(net,clustersWithRoles);

%% plot synaptic weight matrix
if isPlot
    figure(1); clf; 
    %imagesc(S(1:net.N,1:net.N));
    imagesc(S(indexReMap(1:net.N,2),indexReMap(1:net.N,2)));
    xticks([net.indE,net.N]); yticks([net.indE,net.N]);
    xlabel('Weight $(j \rightarrow )$','interpreter','latex','fontsize',22); 
    ylabel('Weight $( \rightarrow i)$','interpreter','latex','fontsize',22);
    title('Weight matrix','interpreter','latex','fontsize',22);
    ax = gca; colorbar(ax,'EastOutside','TickDirection','out'); axis square;
    set(ax,'TickDir','out','box','off');
    drawnow;
end

% eigenvalue calculation takes a long time!
% plot weight matrix eigenvalues
% if isPlot
%     figure(10); clf;
%     eigS = eig(S);
%     plot(real(eigS),imag(eigS),'.b');      
%     xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18); 
%     ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18);
%     title('Weight matrix eigenvalues','interpreter','latex','fontsize',20);
%     set(gca,'TickDir','out','color','none','box','off');
% end

weightSummary = NaN(2*net.Q+1,2*net.Q+1);
for j = 1:size(reMapCell,2)
    for i = 1:size(reMapCell,2)
        weightSummary(i,j) = mean(mean(S(reMapCell{2,i},reMapCell{2,j})));
    end
end

if isPlot
    figure(2); clf;
    imagesc(weightSummary);    
    for i = 1:size(reMapCell,2)
        xline(i-0.5); yline(i-0.5);
        for j = 1:size(reMapCell,2)
            if weightSummary(i,j)<-2.5, textColor = 'w'; else, textColor = 'k'; end
            text(j,i,sprintf('%.3f',weightSummary(i,j)),'Color',textColor,'HorizontalAlignment','center');
        end
    end
    xticks(1:size(reMapCell,2)); yticks(1:size(reMapCell,2));
    xticklabels(reMapCell(1,:)); yticklabels(reMapCell(1,:));
    set(gca,'TickLength',[0,0],'XAxisLocation','top');
    title('Weight matrix: summary (mean weights)','fontsize',22);
    drawnow;
end

%% stimulus
% construct time vector
timev = 0:dt:(totalT+warmup); 
% construct the stimulus
stim.Curve = zeros(1,length(timev));
if strcmp(stim.Profile,'square')     
    stim.Curve(timev>=(warmup+1000*stim.tStart) & timev<=(warmup+1000*stim.tStart+1000*stim.Duration)) = stim.Gain;
elseif strcmp(stim.Profile,'alpha')  
    stim.Curve(timev>=warmup+1000*stim.tStart) = (1/((stim.tauFall*1000)-(stim.tauRise*1000)))*(exp(-(0:dt:(totalT-1000*stim.tStart))/(stim.tauFall*1000))-exp(-(0:dt:(totalT-1000*stim.tStart))/(stim.tauRise*1000)));
    stim.Curve = stim.Curve/max(abs(stim.Curve))*stim.Gain;
end
% get indices of stimulus-selective neurons
stim.sucInd = getStimInd(clustSuc,overlappingIndicesE,stim.frac); % index of sucrose cluster neurons that respond to the stimulus
stim.quiInd = getStimInd(clustQui,overlappingIndicesE,stim.frac); % index of quinine cluster neurons that respond to the stimulus
stim.malInd = getStimInd(clustMal,overlappingIndicesE,stim.frac); % index of maltose cluster neurons that respond to the stimulus
stim.octInd = getStimInd(clustOct,overlappingIndicesE,stim.frac); % index of octaacetate cluster neurons that respond to the stimulus
stim.lInd = sort([stim.sucInd,stim.quiInd]);
stim.rInd = sort([stim.malInd,stim.octInd]);

if isPlot
    figure(3); clf;
    plot((timev-warmup)/1000,stim.Curve*100,'linewidth',2.5,'color','k');
    %plot((0:dt:totalT)/1000,stim.Curve(warmup/dt+1:end)*100,'linewidth',2.5,'color','k');
    xlabel('Time [s]'); ylabel('Stimulus Gain [%]'); title({'Stimulus Profile';''});
    set(gca,'fontsize',10,'TickDir','out','Color','none','box','off');
    drawnow;
end

%% action gate (intended to mimic moving up of lateral spouts)
if isActionGate, gate.Switch = 'on'; else, gate.Switch = 'off'; end
gate.StartTimes = gate.StartMean + gate.Noise*(2*rand(numTrials,1)-1);

if isPlot && isActionGate
    figure(4); clf;
    curve = getActionGate(gate.StartMean,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    plot((timev-warmup)/1000,curve,'linewidth',2.5,'color','k'); hold on;
    %plot((0:dt:totalT)/1000,curve(warmup/dt+1:end),'linewidth',2.5,'color','k'); hold on;
    curve = getActionGate(gate.StartMean-gate.Noise,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    plot((timev-warmup)/1000,curve,'linewidth',1.5,'color','k','linestyle',':'); hold on;
    %plot((0:dt:totalT)/1000,curve(warmup/dt+1:end),'linewidth',1.5,'color','k','linestyle',':'); hold on;
    curve = getActionGate(gate.StartMean+gate.Noise,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    plot((timev-warmup)/1000,curve,'linewidth',1.5,'color','k','linestyle',':'); hold on;
    %plot((0:dt:totalT)/1000,curve(warmup/dt+1:end),'linewidth',1.5,'color','k','linestyle',':'); hold on;
    xlabel('Time [s]'); ylabel('Action Gate'); title({'Action Gate Profile';''});
    ylim([0,1]);
    set(gca,'fontsize',10,'TickDir','out','Color','none','box','off');
    drawnow;
end

%% simulated optogenetic silencing (temporally-controlled increased input to I neurons)
if isSilencing
    silence.curve = getSilencing(silence.start,silence.duration,warmup,totalT,dt);
else
    silence.curve = '';
end

if isPlot && isSilencing
    figure(11); clf;
    plot((timev-warmup)/1000,100*(silence.curve-1),'linewidth',2.5,'color','k');
    %plot((0:dt:totalT)/1000,100*(silence.curve(warmup/dt+1:end)-1),'linewidth',2.5,'color','k');
    xlabel('Time [s]'); ylabel('% Increase Input to I'); title({'Optogenetic Silencing';''});
    set(gca,'fontsize',10','TickDir','out','color','none','box','off');
    drawnow;
end

%% main simulations: baseline
timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nSimulating baseline trial...');

gate.Curve = getActionGate(gate.StartMean,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);

[spkdata,~] = simulation(net,neu,stim,inds,'',gate.Switch,gate.Curve,S,isModifyWeights,warmup,BinW,timev,'',false,'');

fprintf('Done. ');
toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract spike data
firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end

% raster plot: baseline
if isPlot
    figure(5); 
    clf; hold all;
    rasterPlot(firings,inds,indexReMap_inv,clustersWithRoles,totalT,net);
    xlim(([timev(1),timev(end)]-warmup)/1000);
    title('Baseline trial (no stimulus)');
    drawnow;
end

%% main simulations: trials
firings_all = cell(numTrials,1);
inputs_all = cell(numTrials,1);
firingRates = cell(numTrials,6);
stimuli = [repmat({'Sucrose'},1,ceil(numTrials/4)),... 
           repmat({'Maltose'},1,ceil(numTrials/4)),...
           repmat({'Quinine'},1,ceil(numTrials/4)),...
           repmat({'Octaacetate'},1,ceil(numTrials/4))];
       
for iter = 1:numTrials

    timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nSimulating trial %i...',iter);
    gate.Curve = getActionGate(gate.StartTimes(iter),gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    [spkdata,currdata] = simulation(net,neu,stim,inds,stimuli{iter},gate.Switch,gate.Curve,S,isModifyWeights,...
        warmup,BinW,timev,silence.curve,saveInputs,'');

    fprintf('Done. ');
    toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % extract spikes and calculate firing rates
    firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
    firings_all{iter} = firings;
    bins = -warmup/1000:binSize:totalT/1000;
    [~,frLStim] = spikes2FR(firings,bins,[inds.sucE,inds.quiE]); firingRates{iter,1} = frLStim;
    [~,frRStim] = spikes2FR(firings,bins,[inds.malE,inds.octE]); firingRates{iter,2} = frRStim;
    [~,frLCue] = spikes2FR(firings,bins,inds.cueLE); firingRates{iter,3} = frLCue;
    [~,frRCue] = spikes2FR(firings,bins,inds.cueRE); firingRates{iter,4} = frRCue;
    [~,frLAct] = spikes2FR(firings,bins,inds.actLE); firingRates{iter,5} = frLAct;
    [~,frRAct] = spikes2FR(firings,bins,inds.actRE); firingRates{iter,6} = frRAct;
    if saveInputs, inputs_all{iter} = currdata; end

    % raster plot: trial
    if isPlot
        % new figure for each trial if there won't be too many
        if numTrials<=4, figure(iter+5); else, figure(6); end 
        clf; hold all;
        rasterPlot(firings,inds,indexReMap_inv,clustersWithRoles,totalT,net);
        xlim(([timev(1),timev(end)]-warmup)/1000);
        title(sprintf('Trial: %i, Stimulus: %s',iter,stimuli{iter}));
        drawnow;
    end

end

%% threshold firing rates and score performance
corrChoices = reshape(double(cellfun(@(x)ismember(x,{'Maltose','Octaacetate'}),stimuli)),[],1);
firingRates_thresh = cell(numTrials,6);
frLAct_max = max(cellfun(@max,firingRates(:,5)));
frRAct_max = max(cellfun(@max,firingRates(:,6)));
frLCue_max = max(cellfun(@max,firingRates(:,3)));
frRCue_max = max(cellfun(@max,firingRates(:,4)));
for TR = 1:numTrials   
    firingRates_thresh{TR,5} = threshold(firingRates{TR,5},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
    firingRates_thresh{TR,6} = threshold(firingRates{TR,6},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
    firingRates_thresh{TR,3} = threshold(firingRates{TR,3},stateThresh*max(frLCue_max,frRCue_max),stateMinTime,binSize);
    firingRates_thresh{TR,4} = threshold(firingRates{TR,4},stateThresh*max(frLCue_max,frRCue_max),stateMinTime,binSize);
end

choices = NaN(numTrials,1);
bin1 = floor((warmup/1000+stim.tStart+bound1)/binSize)+1; 
bin2 = ceil((warmup/1000+stim.tStart+bound2)/binSize);
for trial = 1:numTrials
    if any(firingRates_thresh{trial,5}(bin1:bin2)) && ~any(firingRates_thresh{trial,6}(bin1:bin2))
        choices(trial) = 0;
    elseif ~any(firingRates_thresh{trial,5}(bin1:bin2)) && any(firingRates_thresh{trial,6}(bin1:bin2))
        choices(trial) = 1;
    end
end

fprintf('\nPerformance by action cluster activation from %.1f to %.1f s (correct/performed/total): %i/%i/%i\n\n',...
        bound1,bound2,sum(choices==corrChoices),sum(~isnan(choices)),numTrials);

toc(timerStart_total);

%% local function definitions
function P_MATRIX = getP(NET,ALL_INDICES_E,ALL_INDICES_I)
% remove connections from synaptic weight matrix according to connectivity
%     by constructing P_MATRIX
% a value of 0 in P_MATRIX marks the corresponding connection in the
%     synaptic weight matrix for removal
% this is a 'fixed random' connectivity scheme, where connections are
%     random but there is a fixed total number of connections from each
%     presynaptic entity to each postsynaptic neuron
% importantly, it takes overlaps into account, meaning that a presynaptic
%     cluster with overlaps is treated as two presynaptic entities, the
%     portion overlapping some other cluster and the portion not in overlap
P_MATRIX = zeros(NET.N,NET.N);
clusters1 = arrayfun(@(x)x.cluster,ALL_INDICES_E);
clusters2 = arrayfun(@(x)x.cluster,ALL_INDICES_I);
if length(clusters1)~=length(clusters2) || length(intersect(clusters1,clusters2))~=length(clusters1)
    error('Arguments for ''ALL_INDICES_E'' and ''ALL_INDICES_I'' must refer to the same set of clusters.');
end
overlappingClusters = clusters1;
for i = 1:NET.Ne
    % clusters->E
    for j = 1:NET.Q
        if ~ismember(j,overlappingClusters)
            % Ec->E
            idx = NET.indE(j):(NET.indE(j+1)-1); % presynaptic indices
            indexLength = length(idx);
            idx = idx(randperm(indexLength)); % scramble
            idx(idx==i) = []; % remove possible self-connections
            P_MATRIX(i,idx(1:round(NET.Pee*indexLength))) = 1; % retain appropriate amount
            % Ic->E
            idx = NET.indI(j):NET.indI(j+1)-1;
            indexLength = length(idx);
            idx = idx(randperm(indexLength));
            P_MATRIX(i,idx(1:round(NET.Pei*indexLength))) = 1; 
        else
            % Ec->E (same as above, but done independently for overlapping/nonoverlapping subclusters)
            allIndices = ALL_INDICES_E(arrayfun(@(x)x.cluster==j,ALL_INDICES_E));
            numConns_overlap = round(NET.Pee*length(allIndices.overlap));
            numConns_nonoverlap = round(NET.Pee*length(allIndices.all))-numConns_overlap;
            % Ec(overlap)->E
            idx1 = allIndices.overlap;
            indexLength1 = length(idx1);
            idx1 = idx1(randperm(indexLength1));
            idx1(idx1==i) = [];
            P_MATRIX(i,idx1(1:numConns_overlap)) = 1;
            % Ec(nonoverlap)->E
            idx2 = allIndices.nonoverlap;
            indexLength2 = length(idx2);
            idx2 = idx2(randperm(indexLength2));
            idx2(idx2==i) = [];
            P_MATRIX(i,idx2(1:numConns_nonoverlap)) = 1; 
            % Ic->E (same as above, but done independently for overlapping/nonoverlapping subclusters)
            allIndices = ALL_INDICES_I(arrayfun(@(x)x.cluster==j,ALL_INDICES_I));
            numConns_overlap = round(NET.Pei*length(allIndices.overlap));
            numConns_nonoverlap = round(NET.Pei*length(allIndices.all))-numConns_overlap;
            % Ic(overlap)->E
            idx1 = allIndices.overlap;
            indexLength1 = length(idx1);
            idx1 = idx1(randperm(indexLength1));
            P_MATRIX(i,idx1(1:numConns_overlap)) = 1; 
            % Ic(nonoverlap)->E
            idx2 = allIndices.nonoverlap;
            indexLength2 = length(idx2);
            idx2 = idx2(randperm(indexLength2));
            P_MATRIX(i,idx2(1:numConns_nonoverlap)) = 1;
        end
    end
    % Eb->E 
    idx = NET.indE(NET.Q+1):NET.Ne;
    indexLength = length(idx);
    idx = idx(randperm(indexLength));
    idx(idx==i) = [];
    P_MATRIX(i,idx(1:round(NET.Pee*indexLength))) = 1; 
end
for i = 1:NET.Ni
    % clusters->I
    for j = 1:NET.Q
        if ~ismember(j,overlappingClusters)
            % Ec->I
            idx = NET.indE(j):(NET.indE(j+1)-1);
            indexLength = length(idx);
            idx = idx(randperm(indexLength));
            P_MATRIX(i+NET.Ne,idx(1:round(NET.Pie*indexLength))) = 1; 
            % Ic->I
            idx = NET.indI(j):NET.indI(j+1)-1;
            indexLength = length(idx);
            idx = idx(randperm(indexLength));
            idx(idx==i+NET.Ne) = [];
            P_MATRIX(i+NET.Ne,idx(1:round(NET.Pii*indexLength))) = 1; 
        else
            % Ec->I
            allIndices = ALL_INDICES_E(arrayfun(@(x)x.cluster==j,ALL_INDICES_E));
            numConns_overlap = round(NET.Pie*length(allIndices.overlap));
            numConns_nonoverlap = round(NET.Pie*length(allIndices.all))-numConns_overlap;
            % Ec(overlap)->I
            idx1 = allIndices.overlap;
            indexLength1 = length(idx1);
            idx1 = idx1(randperm(indexLength1));
            P_MATRIX(i+NET.Ne,idx1(1:numConns_overlap)) = 1; 
            % Ec(nonoverlap)->I
            idx2 = allIndices.nonoverlap;
            indexLength2 = length(idx2);
            idx2 = idx2(randperm(indexLength2));
            P_MATRIX(i+NET.Ne,idx2(1:numConns_nonoverlap)) = 1; 
            % Ic->I
            allIndices = ALL_INDICES_I(arrayfun(@(x)x.cluster==j,ALL_INDICES_I));
            numConns_overlap = round(NET.Pii*length(allIndices.overlap));
            numConns_nonoverlap = round(NET.Pii*length(allIndices.all))-numConns_overlap;
            % Ic(overlap)->I
            idx1 = allIndices.overlap;
            indexLength1 = length(idx1);
            idx1 = idx1(randperm(indexLength1));
            idx1(idx1==i+NET.Ne) = [];
            P_MATRIX(i+NET.Ne,idx1(1:numConns_overlap)) = 1; 
            % Ic(nonoverlap)->I
            idx2 = allIndices.nonoverlap;
            indexLength2 = length(idx2);
            idx2 = idx2(randperm(indexLength2));
            idx2(idx2==i+NET.Ne) = [];
            P_MATRIX(i+NET.Ne,idx2(1:numConns_nonoverlap)) = 1; 
        end
    end
    % Eb->I 
    idx = NET.indE(NET.Q+1):NET.Ne;
    indexLength = length(idx);
    idx = idx(randperm(indexLength));
    P_MATRIX(i+NET.Ne,idx(1:round(NET.Pie*indexLength))) = 1; 
end
end

function MATRIX_NEW = overlap(CLUST1,CLUST2,FACTOR,MATRIX,INDS_E,INDS_I,NET,J_ME,J_PE,J_MI,J_PI)
% take a MATRIX and 'overlap' CLUST1 and CLUST2
% all indexing information for the clusters, including which neurons are
%     and are not supposed to be in the overlapping subcluster, is
%     contained in INDS_E and INDS_I
% the changed weights are determined by FACTOR, NET, J_ME, J_PE, J_MI, and
%     J_PI
J_PPE = J_PE;
J_PE = (1-FACTOR)*J_ME + FACTOR*J_PPE;
J_PPI = J_PI;
J_PI = (1-FACTOR)*J_MI + FACTOR*J_PPI;
J_EE = NET.Jee;
J_EI = NET.Jei;
J_IE = NET.Jie;
J_II = NET.Jii;
S = MATRIX;
for i = [CLUST1,CLUST2]
    for j = [CLUST1,CLUST2]
        clustIndEPost = arrayfun(@(x)x.cluster==i,INDS_E);
        clustIndIPost = arrayfun(@(x)x.cluster==i,INDS_I);
        clustIndEPre = arrayfun(@(x)x.cluster==j,INDS_E);
        clustIndIPre = arrayfun(@(x)x.cluster==j,INDS_I);
        blockDefaultEE = S(INDS_E(clustIndEPost).all,INDS_E(clustIndEPre).all);
        blockDefaultEI = S(INDS_E(clustIndEPost).all,INDS_I(clustIndIPre).all);
        blockDefaultIE = S(INDS_I(clustIndIPost).all,INDS_E(clustIndEPre).all);
        blockDefaultII = S(INDS_I(clustIndIPost).all,INDS_I(clustIndIPre).all);
        if i==j
            blockDefaultEE(blockDefaultEE~=0) = J_PE*J_EE;
            blockDefaultEI(blockDefaultEI~=0) = J_PI*J_EI;
            blockDefaultIE(blockDefaultIE~=0) = J_PI*J_IE;
            blockDefaultII(blockDefaultII~=0) = J_PI*J_II;
        else
            blockDefaultEE(blockDefaultEE~=0) = J_ME*J_EE;
            blockDefaultEI(blockDefaultEI~=0) = J_MI*J_EI;
            blockDefaultIE(blockDefaultIE~=0) = J_MI*J_IE;
            blockDefaultII(blockDefaultII~=0) = J_MI*J_II;
        end
        S(INDS_E(clustIndEPost).all,INDS_E(clustIndEPre).all) = blockDefaultEE;
        S(INDS_E(clustIndEPost).all,INDS_I(clustIndIPre).all) = blockDefaultEI;
        S(INDS_I(clustIndIPost).all,INDS_E(clustIndEPre).all) = blockDefaultIE;
        S(INDS_I(clustIndIPost).all,INDS_I(clustIndIPre).all) = blockDefaultII;
        blockOverlapEE = S(INDS_E(clustIndEPost).overlap,INDS_E(clustIndEPre).overlap);
        blockOverlapEI = S(INDS_E(clustIndEPost).overlap,INDS_I(clustIndIPre).overlap);
        blockOverlapIE = S(INDS_I(clustIndIPost).overlap,INDS_E(clustIndEPre).overlap);
        blockOverlapII = S(INDS_I(clustIndIPost).overlap,INDS_I(clustIndIPre).overlap);
        blockOverlapEE(blockOverlapEE~=0) = J_PPE*J_EE;
        blockOverlapEI(blockOverlapEI~=0) = J_PPI*J_EI;
        blockOverlapIE(blockOverlapIE~=0) = J_PPI*J_IE;
        blockOverlapII(blockOverlapII~=0) = J_PPI*J_II;
        S(INDS_E(clustIndEPost).overlap,INDS_E(clustIndEPre).overlap) = blockOverlapEE;
        S(INDS_E(clustIndEPost).overlap,INDS_I(clustIndIPre).overlap) = blockOverlapEI;
        S(INDS_I(clustIndIPost).overlap,INDS_E(clustIndEPre).overlap) = blockOverlapIE;
        S(INDS_I(clustIndEPost).overlap,INDS_I(clustIndIPre).overlap) = blockOverlapII;
    end
end
MATRIX_NEW = S;
end

function MATRIX_NEW = modifyWeights(POST_IND,PRE_IND,FACTOR,PROB,MATRIX)
% take a block of MATRIX, defined by POST_IND and PRE_IND, and multiply a
%     fraction (PROB) of the entries by FACTOR
% done row by row to ensure uniformity (i.e. all postsynaptic neurons
%     receive the same number of modified inputs)
S = MATRIX;
block = S(POST_IND,PRE_IND); 
for i = 1:size(block,1)
    idx = find(block(i,:)~=0); 
    idx = idx(randperm(length(idx))); 
    idx = sort(idx(1:round(PROB*length(idx))));
    block(i,idx) = FACTOR*block(i,idx);
end
S(POST_IND,PRE_IND) = block;
MATRIX_NEW = S;
end

function [REMAP_CELL,INDEX_REMAP,INDEX_REMAP_INV] = reMapIndex(NET,CLUSTS)
% create mappings between neuron index in terms of the network and neuron
%     index in terms of how you want to plot them
% the function assumes you want to plot neurons in the following order:
%     sucrose E, quinine E, maltose E, octaacetate E, E clusters without
%     roles, background E neurons, sucrose I, quinine I, maltose I,
%     octaacetate I, I clusters without roles
% REMAP_CELL saves each neuron group's network indices in row 2 and plot
%     indices in row 3 (row 1 is the group label)
% INDEX_REMAP has an ordered list of plot indices in its first column and
%     the corresponding network indices in the second column
% INDEX_REMAP_INV has an ordered list of network indices in its first
%     column and the corresponding plot indices in the second column
% so you can easily convert between the two indices like:
%     index_plot = INDEX_REMAP_INV(index_network,2);
%     index_network = INDEX_REMAP(index_plot,2);
INDEX_REMAP = zeros(NET.N,2);
INDEX_REMAP(:,1) = (1:(NET.N))';
count = 1; count2 = 0;
REMAP_CELL = cell(3,2*NET.Q+1);
% sucrose E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(1)):(NET.indE(CLUSTS(1)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Suc_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% quinine E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(2)):(NET.indE(CLUSTS(2)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Qui_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% maltose E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(3)):(NET.indE(CLUSTS(3)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Mal_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% octaacetate E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(4)):(NET.indE(CLUSTS(4)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Oct_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% cue L E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(5)):(NET.indE(CLUSTS(5)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'CueL_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% cue R E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(6)):(NET.indE(CLUSTS(6)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'CueR_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% action L E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(7)):(NET.indE(CLUSTS(7)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'ActL_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% action R E
count2 = count2 + 1;
ind1 = (NET.indE(CLUSTS(8)):(NET.indE(CLUSTS(8)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'ActR_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% other E clusters
for clust = setdiff(1:NET.Q,CLUSTS)
    count2 = count2 + 1;
    ind1 = (NET.indE(clust):(NET.indE(clust+1)-1))';
    ind2 = (count:(count+length(ind1)-1))';
    REMAP_CELL{1,count2} = 'Other_e';
    REMAP_CELL{2,count2} = ind1;
    REMAP_CELL{3,count2} = ind2;
    count = ind2(end)+1;
    INDEX_REMAP(ind2,2) = ind1;
end
% background E neurons
count2 = count2 + 1;
ind1 = (NET.indE(NET.Q+1):(NET.indE(NET.Q+2)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'B_e';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% sucrose I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(1)):(NET.indI(CLUSTS(1)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Suc_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% quinine I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(2)):(NET.indI(CLUSTS(2)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Qui_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% maltose I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(3)):(NET.indI(CLUSTS(3)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Mal_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% octaacetate I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(4)):(NET.indI(CLUSTS(4)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'Oct_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% cue L I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(5)):(NET.indI(CLUSTS(5)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'CueL_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% cue R I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(6)):(NET.indI(CLUSTS(6)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'CueR_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% action L I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(7)):(NET.indI(CLUSTS(7)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'ActL_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% action R I
count2 = count2 + 1;
ind1 = (NET.indI(CLUSTS(8)):(NET.indI(CLUSTS(8)+1)-1))';
ind2 = (count:(count+length(ind1)-1))';
REMAP_CELL{1,count2} = 'ActR_i';
REMAP_CELL{2,count2} = ind1;
REMAP_CELL{3,count2} = ind2;
count = ind2(end)+1;
INDEX_REMAP(ind2,2) = ind1;
% other I clusters
for clust = setdiff(1:NET.Q,CLUSTS)
    count2 = count2 + 1;
    ind1 = (NET.indI(clust):(NET.indI(clust+1)-1))';
    ind2 = (count:(count+length(ind1)-1))';
    REMAP_CELL{1,count2} = 'Other_i';
    REMAP_CELL{2,count2} = ind1;
    REMAP_CELL{3,count2} = ind2;
    count = ind2(end)+1;
    INDEX_REMAP(ind2,2) = ind1;
end
[~,temp] = sort(INDEX_REMAP(:,2));
INDEX_REMAP_INV = [INDEX_REMAP(temp,2),INDEX_REMAP(temp,1)];
end

function INDEX = getStimInd(CLUSTER,ALL_INDICES,FRACTION)
% get the INDEX of stimulus-responsive neurons in a stimulus-responsive
%     CLUSTER
% randomly select a FRACTION of neurons from the overlapping and
%     nonoverlapping pools, specified in ALL_INDICES
INDEX = [];
temp = arrayfun(@(x)x.cluster==CLUSTER,ALL_INDICES);
pool = ALL_INDICES(temp).overlap;
INDEX = [INDEX, pool(randperm(length(pool),round(FRACTION*length(pool))))];
pool = ALL_INDICES(temp).nonoverlap;
INDEX = [INDEX, pool(randperm(length(pool),round(FRACTION*length(pool))))];
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

function SILENCING_CURVE = getSilencing(START,DURATION,STRENGTH,WARMUP,TOTAL_T,DT,varargin)
% Construct simulated optogenetic silencing time course as a square pulse
% input
% Note: if 'center' is given as the last argument to the function, the
% first argument will be used as the center time of the pulse, rather than
% as the start time of the pulse (which is the default behavior ... this is
% what happens if you leave the last argument out)

timeV = 0:DT:(WARMUP+TOTAL_T);
if isempty(varargin)
    % default behavior: use first argument as silencing start time
    start = START;
elseif strcmp(varargin{1},'center')
    % optional behavior: use first argument as silencing center time
    start = START - DURATION/2;
else
    error('Improper syntax. The last argument can only be ''center'' or omitted.');
end
SILENCING_CURVE = ones(1,length(timeV));
SILENCING_CURVE((timeV>=WARMUP+1000*start)&(timeV<=WARMUP+1000*start+1000*DURATION)) = 1 + STRENGTH/100;
end

function [SPIKE_DATA,INPUT_DATA] = simulation(NET,NEU,STIM,INDS,STIMULUS,GATE_SWITCH,GATE_CURVE,...
    S,IS_MODIFY_WEIGHTS,WARMUP,BIN_W,TIME_V,SILENCING,SAVE_INPUTS,V0)
% initialization
q = 1; fired = []; binfired = [];
if isempty(V0)
    V = 0 + 4*randn(NET.N,1); % initial membrane potential
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

function rasterPlot(SPIKE_DATA,INDS,INDEX_REMAP_INV,CLUSTERS_WITH_ROLES,TOTAL_T,NET)
% plot E taste clusters
ind = [INDS.sucE,INDS.quiE,INDS.malE,INDS.octE];
temp = ismember(SPIKE_DATA(:,2),ind);
plot(SPIKE_DATA(temp,1)/1000,INDEX_REMAP_INV(SPIKE_DATA(temp,2),2),'k.','markersize',0.1);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.sucE,2)),'Suc','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.quiE,2)),'Qui','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.malE,2)),'Mal','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.octE,2)),'Oct','VerticalAlignment','middle','fontsize',12);
% plot E cue left and action left clusters
ind = [INDS.cueLE,INDS.actLE];
temp = ismember(SPIKE_DATA(:,2),ind);
plot(SPIKE_DATA(temp,1)/1000,INDEX_REMAP_INV(SPIKE_DATA(temp,2),2),'g.','markersize',0.1);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.cueLE,2)),'Cue L','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.actLE,2)),'Act L','VerticalAlignment','middle','fontsize',12);
% plot E cue right and action right clusters
ind = [INDS.cueRE,INDS.actRE];
temp = ismember(SPIKE_DATA(:,2),ind);
plot(SPIKE_DATA(temp,1)/1000,INDEX_REMAP_INV(SPIKE_DATA(temp,2),2),'b.','markersize',0.1);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.cueRE,2)),'Cue R','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,mean(INDEX_REMAP_INV(INDS.actRE,2)),'Act R','VerticalAlignment','middle','fontsize',12);
% plot remaining E clusters
for clust = setdiff(1:NET.Q,CLUSTERS_WITH_ROLES)
    ind = NET.indE(clust):(NET.indE(clust+1)-1);
    temp = ismember(SPIKE_DATA(:,2),ind);
    plot(SPIKE_DATA(temp,1)/1000,INDEX_REMAP_INV(SPIKE_DATA(temp,2),2),'k.','markersize',0.1);
end
% plot background E neurons
ind = NET.indE(NET.Q+1):(NET.indE(NET.Q+2)-1);
temp = ismember(SPIKE_DATA(:,2),ind);
plot(SPIKE_DATA(temp,1)/1000,INDEX_REMAP_INV(SPIKE_DATA(temp,2),2),'k.','markersize',0.1);
% plot I neurons
ind = (NET.Ne+1):NET.N;
temp = ismember(SPIKE_DATA(:,2),ind);
plot(SPIKE_DATA(temp,1)/1000,INDEX_REMAP_INV(SPIKE_DATA(temp,2),2),'r.','markersize',0.1);
% visual formatting
yline(NET.Ne+0.5,'r');
yline(NET.N+0.5,'r');
set(gca,'fontsize',18,'tickdir','out','color','none','box','off');
yticks([]);
xlim([0, TOTAL_T/1000]); ylim([0, NET.N+1]);
end

function [T_VECTOR,FR_VECTOR] = spikes2FR(SPIKE_DATA,BIN_EDGES,INDS)
binSize = diff(BIN_EDGES(1:2));
spikes = SPIKE_DATA(ismember(SPIKE_DATA(:,2),INDS),:);
FR_VECTOR = zeros(1,length(BIN_EDGES)-1);
for i = 1:length(BIN_EDGES)-1
    FR_VECTOR(i) = sum(spikes(:,1)/1000>=BIN_EDGES(i)&spikes(:,1)/1000<BIN_EDGES(i+1))/length(INDS)/binSize;
end
T_VECTOR = BIN_EDGES(1:end-1)+binSize/2;
end

function CURVE_BINARIZED = threshold(CURVE,THRESHOLD_VALUE,MIN_TIME,BIN_SIZE)
curve_binarized = CURVE>THRESHOLD_VALUE;
onsets = find(diff([0 curve_binarized 0])==1); 
offsets = find(diff([0 curve_binarized 0])==-1);
inds = find(BIN_SIZE*(offsets-onsets)<MIN_TIME);
if ~isempty(inds)
    for ind = inds
        curve_binarized(onsets(ind):min(offsets(ind)-1,length(curve_binarized))) = 0;
    end
end
CURVE_BINARIZED = curve_binarized;
end
