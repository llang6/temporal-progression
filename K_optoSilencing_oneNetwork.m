%% optoSilencing_oneNetwork
%
% This script loads a network and simulates 100 trials with a specified
% silencing configuration. 
%
% We used network 8/10 (a.k.a. 252) since its performance in the absence of
% silencing was closest to average (Figure 6b, Control).
% 
% All spiking data is saved, and then an HMM is fitted.
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Simulation/V0.mat
% /RawData/Simulation/simulationDataX_strongStim.mat for X in 245:254
% /RawData/Simulation/simulationDataX_noSilencing_0_sameV0.mat for X in 245:254 (created internally)
% /RawData/Simulation/simulationData252_noSilencing_0_sameV0.mat (created internally)
% /RawData/Simulation/simulationData252_beginningSilencing_25_sameV0.mat (created internally)
% /RawData/Simulation/simulationData252_beginningSilencing_100_sameV0.mat (created internally)
% /RawData/Simulation/simulationData252_cueOnsetSilencing_25_sameV0.mat (created internally)
% /RawData/Simulation/simulationData252_cueOnsetSilencing_100_sameV0.mat (created internally)
% /RawData/Simulation/simulationData252_middleSilencing_25_sameV0.mat (created internally)
% /RawData/Simulation/simulationData252_middleSilencing_100_sameV0.mat (created internally)
% /ProcessedData/Simulation/spikes_sim252_allConditions.mat (created internally)
% /ProcessedData/Simulation/spikes_sim252_allConditions_multispikesremoved.mat (created internally)
% /ProcessedData/Simulation/win_train_sim252_allConditions.mat (created internally)
%
% requires access to functions:
% loadVar (in +fun)
% getSilencing (in +fun)
% getActionGate (in +fun)
% simulation (in +fun)
% rasterPlot (in +fun)
% spikes2FR (in +fun)
% threshold (in +fun)
% funHMM (in +hmm; requires several other functions in +hmm and +aux)

%% parameters
network_session = 8;
silence_at = 'cueOnset'; % 'no','beginning','cueOnset','middle'
silence_str = 100; % 0, 25, 100
isPlot = true;

%% simulations
% load initial conditions and parameters
V0 = fun.loadVar(sprintf('%s/RawData/Simulation/V0.mat',pwd));
if strcmp(silence_at,'no')
    data_load = load(sprintf('%s/RawData/Simulation/simulationData%i_strongStim.mat',pwd,244+network_session));
else
    data_load = load(sprintf('%s/RawData/Simulation/simulationData%i_noSilencing_0_sameV0.mat',pwd,244+network_session));
    % get cue cluster onset times from control condition if this is not the
    % control condition simulation
    cueOnsetTimes = data_load.cueOnsetTimes;
end
numTrials = 100;
warmup = 0;
dt = data_load.dt;
totalT = 3000; % trial time, [ms]
timev = 0:dt:(totalT+warmup); 
stim = data_load.stim;
gate = data_load.gate;
net = data_load.net;
neu = data_load.neu;
inds = data_load.inds;
S = data_load.S;
isModifyWeights = data_load.isModifyWeights;
BinW = data_load.BinW; % binning for data-chunking during sim
binSize = data_load.binSize; % binning for firing rate calc
stateThresh = data_load.stateThresh;
stateMinTime = data_load.stateMinTime;
indexReMap_inv = data_load.indexReMap_inv;
clustersWithRoles = data_load.clustersWithRoles;
bound1 = data_load.bound1;
bound2 = data_load.bound2;

% silencing
switch silence_at
    case 'no'
        isSilencing = false;
        silence.start = zeros(numTrials,1);
        silence.duration = 0;
        silence.strength = 0;
    case 'beginning'
        isSilencing = true;
        silence.start = zeros(100,1);
        silence.duration = 0.25;
        silence.strength = silence_str;
    case 'cueOnset'
        isSilencing = true;
        silence.start = cueOnsetTimes-0.125; % relative to stimulus onset, [s]
        silence.duration = 0.25;
        silence.strength = silence_str;
    case 'middle'
        isSilencing = true;
        silence.start = 1.375*ones(100,1);
        silence.duration = 0.25;
        silence.strength = silence_str;
end

if isPlot 
    figure(1); clf;
    subplot(1,2,1);
    plot(timev/1000,100*stim.Curve,'k','linewidth',2.5);
    xlabel('Time [s]'); ylabel('% increased input to E'); 
    title('Stimulus');
    set(gca,'fontsize',10','TickDir','out','color','none','box','off');
    subplot(1,2,2);
    if isSilencing
        temp = fun.getSilencing(mean(silence.start),silence.duration,silence.strength,warmup,totalT,dt);
        plot(timev/1000,100*(temp-1),'k','linewidth',2.5);
        clear temp;
    end
    xlabel('Time [s]'); ylabel('% increased input to I'); 
    title('Silencing');
    set(gca,'fontsize',10','TickDir','out','color','none','box','off');
    drawnow;
end

% setup 
timerStart_total = tic;
rngAtStart = rng;
addpath(pwd);

% main simulations: baseline
timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nSimulating baseline trial...');

gate.Curve = fun.getActionGate(gate.StartMean,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
[spkdata_baseline,~] = fun.simulation(net,neu,stim,inds,'',gate.Switch,gate.Curve,S,isModifyWeights,...
    warmup,BinW,timev,'',false,4*randn(4994,1));

fprintf('Done. ');
toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract spike data
firings = []; for i = 1:length(spkdata_baseline), firings = [firings ; spkdata_baseline{i}]; end

% raster plot: baseline
if isPlot
    figure(2); clf; hold all;
    fun.rasterPlot(firings,inds,indexReMap_inv,clustersWithRoles,totalT,net);
    title('Baseline trial (no stimulus)');
    drawnow;
end

% main simulations: trials
firings_all = cell(numTrials,1);
firingRates = cell(numTrials,8);
stimuli = [repmat({'Sucrose'},1,ceil(numTrials/4)),... 
           repmat({'Maltose'},1,ceil(numTrials/4)),...
           repmat({'Quinine'},1,ceil(numTrials/4)),...
           repmat({'Octaacetate'},1,ceil(numTrials/4))];
       
for iter = 1:numTrials

    timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nSimulating trial %i...',iter);

    gate.Curve = fun.getActionGate(gate.StartTimes(iter),gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    % simulated optogenetic silencing (temporally-controlled increased input to I neurons)
    if isSilencing
        silence.curve = fun.getSilencing(silence.start(iter),silence.duration,silence.strength,warmup,totalT,dt);
    else
        silence.curve = '';
    end
    [spkdata,~] = fun.simulation(net,neu,stim,inds,stimuli{iter},gate.Switch,gate.Curve,S,isModifyWeights,...
        warmup,BinW,timev,silence.curve,false,V0(network_session).V0(:,iter));

    fprintf('Done. ');
    toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % extract spikes and calculate firing rates
    firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
    firings_all{iter} = firings;
    bins = -warmup/1000:binSize:totalT/1000;
    [~,frSuc] = fun.spikes2FR(firings,bins,inds.sucE); firingRates{iter,1} = frSuc;
    [~,frQui] = fun.spikes2FR(firings,bins,inds.quiE); firingRates{iter,2} = frQui;
    [~,frMal] = fun.spikes2FR(firings,bins,inds.malE); firingRates{iter,3} = frMal;
    [~,frOct] = fun.spikes2FR(firings,bins,inds.octE); firingRates{iter,4} = frOct;
    [~,frLCue] = fun.spikes2FR(firings,bins,inds.cueLE); firingRates{iter,5} = frLCue;
    [~,frRCue] = fun.spikes2FR(firings,bins,inds.cueRE); firingRates{iter,6} = frRCue;
    [~,frLAct] = fun.spikes2FR(firings,bins,inds.actLE); firingRates{iter,7} = frLAct;
    [~,frRAct] = fun.spikes2FR(firings,bins,inds.actRE); firingRates{iter,8} = frRAct;
    
    % raster plot: trial
    if isPlot
        figure(3); clf; hold all;
        fun.rasterPlot(firings,inds,indexReMap_inv,clustersWithRoles,totalT,net);
        if silence.duration > 0
            patch([silence.start(iter),silence.start(iter)+silence.duration,...
                silence.start(iter)+silence.duration,silence.start(iter)],...
                [0,0,net.N+1,net.N+1],'y','edgecolor','none','facealpha',0.4);
        end
        title(sprintf('Trial: %i, Stimulus: %s',iter,stimuli{iter}));
        drawnow;
    end

end

% threshold firing rates and score performance
corrChoices = reshape(double(cellfun(@(x)ismember(x,{'Maltose','Octaacetate'}),stimuli)),[],1);
firingRates_thresh = cell(numTrials,8);
frSuc_max = max(cellfun(@max,firingRates(:,1)));
frQui_max = max(cellfun(@max,firingRates(:,2)));
frMal_max = max(cellfun(@max,firingRates(:,3)));
frOct_max = max(cellfun(@max,firingRates(:,4)));
frLCue_max = max(cellfun(@max,firingRates(:,5)));
frRCue_max = max(cellfun(@max,firingRates(:,6)));
frLAct_max = max(cellfun(@max,firingRates(:,7)));
frRAct_max = max(cellfun(@max,firingRates(:,8)));
for TR = 1:numTrials    
    firingRates_thresh{TR,1} = fun.threshold(firingRates{TR,1},stateThresh*max([frSuc_max,frQui_max,frMal_max,frOct_max]),stateMinTime,binSize);
    firingRates_thresh{TR,2} = fun.threshold(firingRates{TR,2},stateThresh*max([frSuc_max,frQui_max,frMal_max,frOct_max]),stateMinTime,binSize);
    firingRates_thresh{TR,3} = fun.threshold(firingRates{TR,3},stateThresh*max([frSuc_max,frQui_max,frMal_max,frOct_max]),stateMinTime,binSize);
    firingRates_thresh{TR,4} = fun.threshold(firingRates{TR,4},stateThresh*max([frSuc_max,frQui_max,frMal_max,frOct_max]),stateMinTime,binSize);
    firingRates_thresh{TR,5} = fun.threshold(firingRates{TR,5},stateThresh*max([frLCue_max,frRCue_max]),stateMinTime,binSize);
    firingRates_thresh{TR,6} = fun.threshold(firingRates{TR,6},stateThresh*max([frLCue_max,frRCue_max]),stateMinTime,binSize);
    firingRates_thresh{TR,7} = fun.threshold(firingRates{TR,7},stateThresh*max([frLAct_max,frRAct_max]),stateMinTime,binSize);
    firingRates_thresh{TR,8} = fun.threshold(firingRates{TR,8},stateThresh*max([frLAct_max,frRAct_max]),stateMinTime,binSize);
end

choices = NaN(numTrials,1);
bin1 = floor((stim.tStart+bound1)/binSize)+1; 
bin2 = ceil((stim.tStart+bound2)/binSize);
for trial = 1:numTrials
    if any(firingRates_thresh{trial,7}(bin1:bin2)) && ~any(firingRates_thresh{trial,8}(bin1:bin2))
        choices(trial) = 0;
    elseif ~any(firingRates_thresh{trial,7}(bin1:bin2)) && any(firingRates_thresh{trial,8}(bin1:bin2))
        choices(trial) = 1;
    end
end

fprintf('\nStrict performance by action cluster activation from %.1f to %.1f s (correct/performed/total): %i/%i/%i\n\n',...
        binSize*(bin1-1),binSize*bin2,sum(choices==corrChoices),sum(~isnan(choices)),numTrials);

if strcmp(silence_at,'no')
    % get cue cluster onset times from current simulation if this is the
    % control condition simulation
    cueOnsetTimes = 0.5*ones(numTrials,1);
    for i = 1:numTrials
        cueLeftTime = inf;
        cueRightTime = inf;
        if ~isempty(find(firingRates_thresh{i,5},1))
            cueLeftTime = (find(firingRates_thresh{i,5},1)-1)*binSize+binSize/2;
        end
        if ~isempty(find(firingRates_thresh{i,6},1))
            cueRightTime = (find(firingRates_thresh{i,6},1)-1)*binSize+binSize/2;
        end
        if ~all(isinf([cueLeftTime,cueRightTime]))
            cueOnsetTimes(i) = min([cueLeftTime,cueRightTime]);
        else
            fprintf('\nWarning: No cue cluster onset in trial %i. Using default value of 0.5 s\n',i);
        end
    end
end
    
toc(timerStart_total);

% save(sprintf('%s/RawData/Simulation/simulationData%i_%sSilencing_%i_sameV0',pwd,244+network_session,silence_at,silence_str),...
%     '-regexp','^(?!data_load$).');

%% prepare the spiking data for HMM fitting
% this requires simulation files created using the cells above

% parameters
%clearvars; clc;
preStim = 0;
postStim = 3;
nSamp = 1;
nDownsample = 14;

% load data
data_baseline = load(sprintf('%s/RawData/Simulation/simulationData252_noSilencing_0_sameV0.mat',pwd));
data_silenceBeginning_25 = load(sprintf('%s/RawData/Simulation/simulationData252_beginningSilencing_25_sameV0.mat',pwd));
data_silenceBeginning_100 = load(sprintf('%s/RawData/Simulation/simulationData252_beginningSilencing_100_sameV0.mat',pwd));
data_silenceCue_25 = load(sprintf('%s/RawData/Simulation/simulationData252_cueOnsetSilencing_25_sameV0.mat',pwd));
data_silenceCue_100 = load(sprintf('%s/RawData/Simulation/simulationData252_cueOnsetSilencing_100_sameV0.mat',pwd));
data_silenceMiddle_25 = load(sprintf('%s/RawData/Simulation/simulationData252_middleSilencing_25_sameV0.mat',pwd));
data_silenceMiddle_100 = load(sprintf('%s/RawData/Simulation/simulationData252_middleSilencing_100_sameV0.mat',pwd));

% adjust windows
firings_all = data_baseline.firings_all;
ntrials = numel(firings_all);
stimStart = data_baseline.stim.tStart; trialEnd = data_baseline.totalT/1000;
if stimStart - preStim < 0, preStim = stimStart; end
if stimStart + postStim > trialEnd, postStim = trialEnd - stimStart; end
leftBound = (stimStart - preStim)*1000; rightBound = (stimStart + postStim)*1000;

win_train_base = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
win_train_silBeg_25 = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
win_train_silBeg_100 = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
win_train_silCue_25 = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
win_train_silCue_100 = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
win_train_silMid_25 = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
win_train_silMid_100 = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];

% randomly sample neurons (nSamp per E cluster, nDownsample total)
sampNeurons = NaN(nSamp,data_baseline.net.Q); % indices of sampled neurons
nCount = 0;
for cluster = [data_baseline.clustersWithRoles,setdiff(1:data_baseline.net.Q,data_baseline.clustersWithRoles)]
    nCount = nCount + 1;
    indices = (data_baseline.net.indE(cluster)):(data_baseline.net.indE(cluster+1)-1);
    FRs = NaN(length(indices),1);
    for i = 1:length(indices)
        frates = zeros(1,ntrials);
        for j = 1:ntrials
            myArray = firings_all{j};
            mySpikes = myArray(myArray(:,2)==indices(i),1);
            frates(j) = sum(mySpikes>leftBound & mySpikes<rightBound)/(preStim+postStim);
        end
        FRs(i) = mean(frates);
    end
    if sum(FRs>=2)<2 
        fprintf('\nWARNING! Cluster %i has fewer than 2 neurons with firing rates at least 2 Hz. Randomly choosing %i and moving on...\n',cluster,nSamp);  
    else
        indices = indices(FRs>=2);
    end
    rN = indices(randperm(length(indices))); % random sampling
    sampNeurons(:,nCount) = reshape(rN(1:nSamp),[],1);
end
sampNeurons = reshape(sampNeurons,1,[]);
sampNeurons = sampNeurons(sort(randperm(length(sampNeurons),nDownsample))); % random sub-sampling

% get the spike trains
spikes_base = struct();
spikes_silBeg_25 = struct();
spikes_silBeg_100 = struct();
spikes_silCue_25 = struct();
spikes_silCue_100 = struct();
spikes_silMid_25 = struct();
spikes_silMid_100 = struct();

for trial = 1:ntrials
    for neuron = 1:length(sampNeurons)
        f = firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_base(trial,neuron).spk = f(f>=win_train_base(trial,1)&f<=win_train_base(trial,2));
        
        f = data_silenceBeginning_25.firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_silBeg_25(trial,neuron).spk = f(f>=win_train_silBeg_25(trial,1)&f<=win_train_silBeg_25(trial,2));
        
        f = data_silenceBeginning_100.firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_silBeg_100(trial,neuron).spk = f(f>=win_train_silBeg_100(trial,1)&f<=win_train_silBeg_100(trial,2));
        
        f = data_silenceCue_25.firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_silCue_25(trial,neuron).spk = f(f>=win_train_silCue_25(trial,1)&f<=win_train_silCue_25(trial,2));
        
        f = data_silenceCue_100.firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_silCue_100(trial,neuron).spk = f(f>=win_train_silCue_100(trial,1)&f<=win_train_silCue_100(trial,2));
        
        f = data_silenceMiddle_25.firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_silMid_25(trial,neuron).spk = f(f>=win_train_silMid_25(trial,1)&f<=win_train_silMid_25(trial,2));
        
        f = data_silenceMiddle_100.firings_all{trial};
        f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
        spikes_silMid_100(trial,neuron).spk = f(f>=win_train_silMid_100(trial,1)&f<=win_train_silMid_100(trial,2));
    end
end

% concatenate
spikes = [spikes_base;...
    spikes_silBeg_25;...
    spikes_silBeg_100;...
    spikes_silCue_25;...
    spikes_silCue_100;...
    spikes_silMid_25;...
    spikes_silMid_100];

win_train = [win_train_base;...
    win_train_silBeg_25;...
    win_train_silBeg_100;...
    win_train_silCue_25;...
    win_train_silCue_100;...
    win_train_silMid_25;...
    win_train_silMid_100];

% save(sprintf('%s/ProcessedData/Simulation/spikes_sim252_allConditions.mat',pwd),'spikes');
% save(sprintf('%s/ProcessedData/Simulation/win_train_sim252_allConditions.mat',pwd),'win_train');

%% continue preparing data for HMM fitting
% here multiple spikes occurring in the same time bin are removed
% beforehand to guarantee that the algorithm always gets the same result
% when asked to bin the same spike train (normally wouldn't matter but we
% want everything other than silencing to be the same across silencing
% conditions)

clearvars; clc;

spikes = fun.loadVar(sprintf('%s/ProcessedData/Simulation/spikes_sim252_allConditions.mat',pwd));

spikes_temp = struct();
for trial = 1:700
    rng(9995); % set rng before each trial so trials with same initial conditions are also binned the same way
    % code below comes directly from binning functions of fitting algorithm
    % (hmm.Spikes2Seq and hmm.low_bernoulli)
    temp_all = arrayfun(@(x)x.spk,spikes(trial,:),'UniformOutput',false);
    temp_shift = cellfun(@(x)x,temp_all(:),'UniformOutput',false)';
    Bins = ceil(3/0.002);
    gnunits = numel(temp_shift);
    X_temp = zeros(Bins,gnunits);
    % eliminate pre-cue and post-delivery spikes
    a_pos = cellfun(@(x)x(x(:)>0 & x(:)<=3),temp_shift,'UniformOutput',false);
    % turn spike times into ms and ceil them to get indices of Bins
    ind_spk = cellfun(@(x)ceil(x(:)/0.002),a_pos,'UniformOutput',false);
    % set X_temp indices to one for each spike
    for unit = 1:gnunits
        X_temp(ind_spk{unit}(:),unit) = 1;
    end
    % find bins with more than one spike
    ind_many = find(sum(X_temp')>1);
    % # of multiple spikes per bin in bins with more than one spike
    SkipSpikes = sum(X_temp(ind_many,:),2);
    temp_a = a_pos;
    temp_i = ind_spk;
    for b = 1:numel(ind_many)
        u = [];
        % find multiple units firing in that bin
        u = find(X_temp(ind_many(b),:));
        u_ind = randperm(numel(u));
        row = ind_many(b);
        columns = u(u_ind(1:end-1));
        X_temp(row,columns) = 0;
        % remove from original spike train
        for i = reshape(columns,1,[])
            temp1 = temp_a{i};
            temp2 = temp_i{i};
            temp1(temp2==row) = [];
            temp2(temp2==row) = [];
            temp_a{i} = temp1;
            temp_i{i} = temp2;
        end
    end
    for unit = 1:gnunits
        spikes_temp(trial,unit).spk = temp_a{unit};
    end
end
spikes = spikes_temp;
                
% save(sprintf('%s/ProcessedData/Simulation/spikes_sim252_allConditions_multispikesremoved.mat',pwd),'spikes');

%% fit the HMM
spikes = fun.loadVar(sprintf('%s/ProcessedData/Simulation/spikes_sim252_allConditions_multispikesremoved.mat',pwd));
win_train = fun.loadVar(sprintf('%s/ProcessedData/Simulation/win_train_sim252_allConditions.mat',pwd));
[ntrials, gnunits] = size(spikes);
% HMM parameters
MODELSEL = 'BIC'; % 'XVAL', 'AIC'
DATAIN = struct('spikes',spikes,'win',win_train,'METHOD',MODELSEL); 
% pass all to main HMM function
res = hmm.funHMM(DATAIN);
% save file
% save(sprintf('%s/HMMData/Simulation/HMM_sim252_allConditions_multispikesremoved.mat',pwd),'res','spikes','win_train');
