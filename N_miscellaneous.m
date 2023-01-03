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
% /ProcessedData/Simulation/perf_silStr25_silDur250_stimFall160_stimGain200.mat
% /ProcessedData/Simulation/perf_silStr100_silDur250_stimFall160_stimGain200.mat
%
% requires access to functions:
% loadVar (in +fun)
% distinguishable_colors (in +fun)
% getStateCell (in +fun)
% scoreTrials (in +fun)
% hist2curve (in +fun)
% mySmoothing (in +fun)
% rasterPlotAbbrev (in +fun)
% getActionGate (in +fun)
% simulation (in +fun)

%% setup
addpath(pwd); % gain access to package functions

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
stateData = fun.loadVar(fileName);

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

% How many sessions contain Quality-coding AND Cue-coding states?
numSessions = 0;
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Cue-coding'}))
        numSessions = numSessions + 1;
    end
end
fprintf('\nNumber of sessions with Quality-coding and Cue-coding states: %i\n',numSessions);

% How many sessions contain Quality-coding AND Action-coding states?
numSessions = 0;
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Action-coding'}))
        numSessions = numSessions + 1;
    end
end
fprintf('\nNumber of sessions with Quality-coding and Action-coding states: %i\n',numSessions);

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

%% Additional state counts
% How many sessions/trials contain any Quality-, Cue-, or Action-coding state?
nTrials_total = [];
nTrials_any = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    nTrials_total = [nTrials_total, length(stateData{find(tf1,1),9})];
    tf_q = [false;cellfun(@(x)strcmp(x,'Exclusive Quality-coding'),stateData(2:end,3))];
    trials_q = [stateData{tf1&tf_q,8}];
    tf_c = [false;cellfun(@(x)strcmp(x,'Cue-coding'),stateData(2:end,3))];
    trials_c = [stateData{tf1&tf_c,8}];
    tf_a = [false;cellfun(@(x)strcmp(x,'Action-coding'),stateData(2:end,3))];
    trials_a = [stateData{tf1&tf_c,8}];
    trials_any = unique([trials_q,trials_c,trials_a]);
    nTrials_any = [nTrials_any, length(trials_any)];
end
fprintf('\nOut of %i total trials, %i contain a coding state (Quality, Cue, or Action)\n',sum(nTrials_total),sum(nTrials_any));

% How many contain the following sequences...?

% Quality-coding --> Cue-coding (no Action-coding)
sessionsOfInterest = [];
nTrials_total = [];
nTrials_multistate = [];
nTrials_inorder = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Cue-coding'})) && ...
            ~any(ismember(block,{'Action-coding'}))
        sessionsOfInterest = [sessionsOfInterest, i];
        nTrials_total = [nTrials_total, length(stateData{find(tf1,1),9})];
        tf_q = [false;cellfun(@(x)strcmp(x,'Exclusive Quality-coding'),stateData(2:end,3))];
        trials_q = [stateData{tf1&tf_q,8}];
        tf_c = [false;cellfun(@(x)strcmp(x,'Cue-coding'),stateData(2:end,3))];
        trials_c = [stateData{tf1&tf_c,8}];
        trials_both = intersect(trials_q,trials_c);
        nTrials_multistate = [nTrials_multistate,length(trials_both)];
        inorder_indicator = 0;
        for trial = trials_both
            onsets_q = [stateData{tf1&tf_q,4}];
            onsets_q = onsets_q([stateData{tf1&tf_q,8}]==trial);
            onsets_c = [stateData{tf1&tf_c,4}];
            onsets_c = onsets_c([stateData{tf1&tf_c,8}]==trial);
            if min(onsets_q)<min(onsets_c)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_inorder = [nTrials_inorder,inorder_indicator];    
    end
end
fprintf('\nQuality --> Cue (no Action) progression:');
fprintf('\n    Occurs in %i/17 sessions',length(sessionsOfInterest));
fprintf('\n    How many trials in order? How many multi-state trials? How many trials in total?');
fprintf('\n    %i/%i/%i',nTrials_inorder(1),nTrials_multistate(1),nTrials_total(1));
if length(sessionsOfInterest)>1
    for i = 2:length(sessionsOfInterest)
        fprintf('; %i/%i/%i',nTrials_inorder(i),nTrials_multistate(i),nTrials_total(i));
    end
end
fprintf('\n\n');

% Quality-coding --> Action-coding (no Cue-coding)
sessionsOfInterest = [];
nTrials_total = [];
nTrials_multistate = [];
nTrials_inorder = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Action-coding'})) && ...
            ~any(ismember(block,{'Cue-coding'}))
        sessionsOfInterest = [sessionsOfInterest, i];
        nTrials_total = [nTrials_total, length(stateData{find(tf1,1),9})];
        tf_q = [false;cellfun(@(x)strcmp(x,'Exclusive Quality-coding'),stateData(2:end,3))];
        trials_q = [stateData{tf1&tf_q,8}];
        tf_a = [false;cellfun(@(x)strcmp(x,'Action-coding'),stateData(2:end,3))];
        trials_a = [stateData{tf1&tf_a,8}];
        trials_both = intersect(trials_q,trials_a);
        nTrials_multistate = [nTrials_multistate,length(trials_both)];
        inorder_indicator = 0;
        for trial = trials_both
            onsets_q = [stateData{tf1&tf_q,4}];
            onsets_q = onsets_q([stateData{tf1&tf_q,8}]==trial);
            onsets_a = [stateData{tf1&tf_a,4}];
            onsets_a = onsets_a([stateData{tf1&tf_a,8}]==trial);
            if min(onsets_q)<min(onsets_a)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_inorder = [nTrials_inorder,inorder_indicator];    
    end
end
fprintf('\nQuality --> Action (no Cue) progression:');
fprintf('\n    Occurs in %i/17 sessions',length(sessionsOfInterest));
fprintf('\n    How many trials in order? How many multi-state trials? How many trials in total?');
fprintf('\n    %i/%i/%i',nTrials_inorder(1),nTrials_multistate(1),nTrials_total(1));
if length(sessionsOfInterest)>1
    for i = 2:length(sessionsOfInterest)
        fprintf('; %i/%i/%i',nTrials_inorder(i),nTrials_multistate(i),nTrials_total(i));
    end
end
fprintf('\n\n');

% Cue-coding --> Action-coding (no Quality-coding)
sessionsOfInterest = [];
nTrials_total = [];
nTrials_multistate = [];
nTrials_inorder = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Cue-coding'})) && ...
            any(ismember(block,{'Action-coding'})) && ...
            ~any(ismember(block,{'Exclusive Quality-coding'}))
        sessionsOfInterest = [sessionsOfInterest, i];
        nTrials_total = [nTrials_total, length(stateData{find(tf1,1),9})];
        tf_c = [false;cellfun(@(x)strcmp(x,'Cue-coding'),stateData(2:end,3))];
        trials_c = [stateData{tf1&tf_c,8}];
        tf_a = [false;cellfun(@(x)strcmp(x,'Action-coding'),stateData(2:end,3))];
        trials_a = [stateData{tf1&tf_a,8}];
        trials_both = intersect(trials_c,trials_a);
        nTrials_multistate = [nTrials_multistate,length(trials_both)];
        inorder_indicator = 0;
        for trial = trials_both
            onsets_c = [stateData{tf1&tf_c,4}];
            onsets_c = onsets_c([stateData{tf1&tf_c,8}]==trial);
            onsets_a = [stateData{tf1&tf_a,4}];
            onsets_a = onsets_a([stateData{tf1&tf_a,8}]==trial);
            if min(onsets_c)<min(onsets_a)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_inorder = [nTrials_inorder,inorder_indicator];    
    end
end
fprintf('\nCue --> Action (no Quality) progression:');
fprintf('\n    Occurs in %i/17 sessions',length(sessionsOfInterest));
fprintf('\n    How many trials in order? How many multi-state trials? How many trials in total?');
fprintf('\n    %i/%i/%i',nTrials_inorder(1),nTrials_multistate(1),nTrials_total(1));
if length(sessionsOfInterest)>1
    for i = 2:length(sessionsOfInterest)
        fprintf('; %i/%i/%i',nTrials_inorder(i),nTrials_multistate(i),nTrials_total(i));
    end
end
fprintf('\n\n');

% Quality-coding --> Cue-coding --> Action-coding
sessionsOfInterest = [];
nTrials_total = [];
nTrials_qc = [];
nTrials_qc_inorder = [];
nTrials_qa = [];
nTrials_qa_inorder = [];
nTrials_ca = [];
nTrials_ca_inorder = [];
nTrials_qca = [];
nTrials_qca_inorder = [];
for i = sessions
    tf1 = [false;cellfun(@(x)x==i,stateData(2:end,1))];
    block = stateData(tf1,3);
    if any(ismember(block,{'Exclusive Quality-coding'})) && ...
            any(ismember(block,{'Cue-coding'})) && ...
            any(ismember(block,{'Action-coding'}))
        
        sessionsOfInterest = [sessionsOfInterest, i];
        nTrials_total = [nTrials_total, length(stateData{find(tf1,1),9})];
        tf_q = [false;cellfun(@(x)strcmp(x,'Exclusive Quality-coding'),stateData(2:end,3))];
        trials_q = [stateData{tf1&tf_q,8}];
        tf_c = [false;cellfun(@(x)strcmp(x,'Cue-coding'),stateData(2:end,3))];
        trials_c = [stateData{tf1&tf_c,8}];
        tf_a = [false;cellfun(@(x)strcmp(x,'Action-coding'),stateData(2:end,3))];
        trials_a = [stateData{tf1&tf_a,8}];
        
        trials_qc = setdiff(intersect(trials_q,trials_c),trials_a);
        trials_qa = setdiff(intersect(trials_q,trials_a),trials_c);
        trials_ca = setdiff(intersect(trials_c,trials_a),trials_q);
        trials_qca = intersect(intersect(trials_q,trials_c),trials_a); 
        nTrials_qc = [nTrials_qc,length(trials_qc)];
        nTrials_qa = [nTrials_qa,length(trials_qa)];
        nTrials_ca = [nTrials_ca,length(trials_ca)];
        nTrials_qca = [nTrials_qca,length(trials_qca)];
        
        inorder_indicator = 0;
        for trial = trials_qc
            onsets_q = [stateData{tf1&tf_q,4}];
            onsets_q = onsets_q([stateData{tf1&tf_q,8}]==trial);
            onsets_c = [stateData{tf1&tf_c,4}];
            onsets_c = onsets_c([stateData{tf1&tf_c,8}]==trial);
            if min(onsets_q)<min(onsets_c)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_qc_inorder = [nTrials_qc_inorder,inorder_indicator];  
        
        inorder_indicator = 0;
        for trial = trials_qa
            onsets_q = [stateData{tf1&tf_q,4}];
            onsets_q = onsets_q([stateData{tf1&tf_q,8}]==trial);
            onsets_a = [stateData{tf1&tf_a,4}];
            onsets_a = onsets_a([stateData{tf1&tf_a,8}]==trial);
            if min(onsets_q)<min(onsets_a)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_qa_inorder = [nTrials_qa_inorder,inorder_indicator];
        
        inorder_indicator = 0;
        for trial = trials_ca
            onsets_c = [stateData{tf1&tf_c,4}];
            onsets_c = onsets_c([stateData{tf1&tf_c,8}]==trial);
            onsets_a = [stateData{tf1&tf_a,4}];
            onsets_a = onsets_a([stateData{tf1&tf_a,8}]==trial);
            if min(onsets_c)<min(onsets_a)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_ca_inorder = [nTrials_ca_inorder,inorder_indicator];
        
        inorder_indicator = 0;
        for trial = trials_qca
            onsets_q = [stateData{tf1&tf_q,4}];
            onsets_q = onsets_q([stateData{tf1&tf_q,8}]==trial);
            onsets_c = [stateData{tf1&tf_c,4}];
            onsets_c = onsets_c([stateData{tf1&tf_c,8}]==trial);
            onsets_a = [stateData{tf1&tf_a,4}];
            onsets_a = onsets_a([stateData{tf1&tf_a,8}]==trial);
            if min(onsets_q)<min(onsets_c) && min(onsets_c)<min(onsets_a)
                inorder_indicator = inorder_indicator + 1;
            end
        end
        nTrials_qca_inorder = [nTrials_qca_inorder,inorder_indicator];
        
    end
end
fprintf('\nQuality --> Cue --> Action progression:');
fprintf('\n    Occurs in %i/17 sessions',length(sessionsOfInterest));
fprintf('\n    How many trials in order? How many multi-state trials? How many trials in total?\n');
fprintf('\n    QC: %i/%i/%i',nTrials_qc_inorder(1),nTrials_qc(1),nTrials_total(1));
fprintf('\n    QA: %i/%i/%i',nTrials_qa_inorder(1),nTrials_qa(1),nTrials_total(1));
fprintf('\n    CA: %i/%i/%i',nTrials_ca_inorder(1),nTrials_ca(1),nTrials_total(1));
fprintf('\n    QCA: %i/%i/%i',nTrials_qca_inorder(1),nTrials_qca(1),nTrials_total(1));
if length(sessionsOfInterest)>1
    for i = 2:length(sessionsOfInterest)
        fprintf('\n\n    QC: %i/%i/%i',nTrials_qc_inorder(i),nTrials_qc(i),nTrials_total(i));
        fprintf('\n    QA: %i/%i/%i',nTrials_qa_inorder(i),nTrials_qa(i),nTrials_total(i));
        fprintf('\n    CA: %i/%i/%i',nTrials_ca_inorder(i),nTrials_ca(i),nTrials_total(i));
        fprintf('\n    QCA: %i/%i/%i',nTrials_qca_inorder(i),nTrials_qca(i),nTrials_total(i));
    end
end
fprintf('\n\n');

%% Behavior summary stats
SessionInfo = fun.loadVar(sprintf('%s/RawData/Experiment/SessionInfo.mat',pwd));

IEIs = [];
for i = 1:length(SessionInfo)
    IEIs = [IEIs,SessionInfo(i).TrialEvents.DecisionTimes-SessionInfo(i).TrialEvents.SampleTimes];
end
fprintf('\nIEI mean: %.2f, range: %.2f-%.2f\n',mean(IEIs),min(IEIs),max(IEIs));

num_neurons = [];
frs = cell(21,1);
for i = 1:length(SessionInfo)
    frs{i} = [];
    neu_count = 0;
    sampling_times = SessionInfo(i).TrialEvents.SampleTimes;
    decision_times = SessionInfo(i).TrialEvents.DecisionTimes;
    for j = 1:length(SessionInfo(i).Neuron)
        spike_train = SessionInfo(i).Neuron(j).SpikeTrain;
        spk_count = 0;
        for k = 1:length(sampling_times)
            spk_count = spk_count + sum(spike_train>=(sampling_times(k)-0.1)&spike_train<=(decision_times(k)+0.1));
        end
        fr = spk_count/(sum(diff([sampling_times;decision_times]))+0.2*length(sampling_times));
        frs{i} = [frs{i};fr];
        if  fr >= 2
            neu_count = neu_count + 1;
        end
    end
    num_neurons = [num_neurons,neu_count];
end

sessions_excluded = find(num_neurons<3);
sessions_included = setdiff(1:length(SessionInfo),sessions_excluded);

num_trials = arrayfun(@(x)length(x.SampleTimes),[SessionInfo(:).TrialEvents]);
accuracy = 100*arrayfun(@(x)mean([x.TrialType(:).wasCorrect]),[SessionInfo(:).TrialEvents]);

fprintf('\nNumber of trials (no excluded sessions) mean: %.1f, range: %i-%i\n',mean(num_trials),min(num_trials),max(num_trials));
fprintf('\nNumber of trials (excluded sessions) mean: %.1f, range: %i-%i\n',mean(num_trials(sessions_included)),min(num_trials(sessions_included)),max(num_trials(sessions_included)));

fprintf('\nAccuracy (no excluded sessions) mean: %.1f, range: %.1f-%.1f\n',mean(accuracy),min(accuracy),max(accuracy));
fprintf('\nAccuracy (excluded sessions) mean: %.1f, range: %.1f-%.1f\n',mean(accuracy(sessions_included)),min(accuracy(sessions_included)),max(accuracy(sessions_included)));

fprintf('\nNumber of neurons (no excluded sessions) mean: %.1f, range: %i-%i\n',mean(num_neurons),min(num_neurons),max(num_neurons));
fprintf('\nNumber of neurons (excluded sessions) mean: %.1f, range: %i-%i\n',mean(num_neurons(sessions_included)),min(num_neurons(sessions_included)),max(num_neurons(sessions_included)));

%% HMM state duration distribution
sessions_excluded = [5,7,17,19,20];
sessions = setdiff(1:21,sessions_excluded);
state_durations = [];
for i = 1:length(sessions)
    hmmData = load(sprintf('%s/HMMData/Experiment/HMM_exp%i.mat',pwd,sessions(i)));
    postfit = hmmData.res.hmm_postfit;
    for j = 1:length(postfit)
        seq = postfit(j).sequence;
        if isempty(seq), continue; end
        state_durations = [state_durations,seq(3,:)];
    end
end

figure(1); clf;
histogram(state_durations,'normalization','pdf');
set(gca,'tickdir','out','color','none','box','off');
xlabel('State durations [s]');
ylabel('Probability density');

fprintf('\nMean: %.1f ms\n',mean(state_durations*1000));
fprintf('\nMedian: %.1f ms\n',median(state_durations*1000));
fprintf('\nStdev: %.1f ms\n',sqrt(var(state_durations*1000)));
fprintf('\nRange: %.1f ms --- %.1f ms\n',min(state_durations*1000),max(state_durations*1000));

%% Plot any experimental trial with HMM state overlay (for Figure 1b-d, left)

% Figure 1b (left) is Session 13, Trial 218 (state 3)
% Figure 1c (left) is Session 16, Trial 166 (state 2)
% Figure 1d (left) is Session 2, Trial 115 (state 2)

% ---------------------
session = 2;
trial = 115;
% ---------------------

classifiedData = fun.loadVar(sprintf('%s/HMMData/Experiment/classifiedStates_exp.mat',pwd));
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
colors = hmmData.res.colors;
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

%% Plot classified states within all trials for any experimental session
% ---------------------
session = 12; % 1-21, ignore 5,7,17,19,20
model = ''; % '', '_shuff_circ', '_shuff_swap'
isWarp = 1;
colorByClass = 0;
% ---------------------

classifiedData = fun.loadVar(sprintf('%s/HMMData/Experiment/classifiedStates_exp%s.mat',pwd,model));
hmmData = load(sprintf('%s/HMMData/Experiment/HMM_exp%i%s.mat',pwd,session,model));
win_train = hmmData.win_train;
behaviorData = load(sprintf('%s/RawData/Experiment/SessionInfo.mat',pwd));
SessionInfo = behaviorData.SessionInfo;
trial_stim = {SessionInfo(session).TrialEvents.TrialType(:).TasteID}';
trial_outcome = vertcat(SessionInfo(session).TrialEvents.TrialType(:).wasCorrect);
stimuli_sorted = {'Sucrose','Maltose','Quinine','Octaacetate'};

fig = figure(2); clf; hold all;
cond2 = trial_outcome==1;
counter = 0;
tickMark = 0;
infoMat = NaN(length(cond2),3);
% error trials
for c = 1:length(stimuli_sorted)
    cond1 = cellfun(@(x)strcmp(x,stimuli_sorted{c}),trial_stim);
    trials = find(cond1&~cond2);
    for i = 1:length(trials)
        counter = counter + 1;
        tickMark = tickMark + 1;
        infoMat(counter,:) = [c,tickMark,0];
        seq = hmmData.res.hmm_postfit(trials(i)).sequence;
        IEI = diff(win_train(trials(i),:))-0.2;        
        if isempty(seq), continue; end
        for j = 1:size(seq,2)
            x = [seq(1,j),seq(2,j),seq(2,j),seq(1,j)] + IEI;
            if isWarp
                for k = 1:4
                    if x(k)>=0 && x(k)<=IEI
                       x(k) = x(k)/IEI;
                    elseif x(k)>IEI
                       x(k) = x(k)-IEI+1;
                    end
                end
            end
            y = [tickMark-0.4,tickMark-0.4,tickMark+0.4,tickMark+0.4];
            if colorByClass
                tf = vertcat(classifiedData{2:end,1})==session & vertcat(classifiedData{2:end,2})==seq(4,j);
                if contains(classifiedData{[false;tf],3},'Quality')
                    color = 'r';
                elseif contains(classifiedData{[false;tf],3},'Cue')
                    color = 'c';
                elseif contains(classifiedData{[false;tf],3},'Action')
                    color = 'b';
                else
                    color = [0.7,0.7,0.7];
                end 
            else
                color = hmmData.res.colors(seq(4,j),:);
            end
            patch(x,y,color,'edgecolor','none','facealpha',0.4);
        end
    end
    if ~isempty(trials)
        tickMark = tickMark + 2;
    end
end
tickMark = tickMark + 10;
% correct trials
for c = 1:length(stimuli_sorted)
    cond1 = cellfun(@(x)strcmp(x,stimuli_sorted{c}),trial_stim);
    trials = find(cond1&cond2);
    for i = 1:length(trials)
        counter = counter + 1;
        tickMark = tickMark + 1;
        infoMat(counter,:) = [c,tickMark,1];
        seq = hmmData.res.hmm_postfit(trials(i)).sequence;
        IEI = diff(win_train(trials(i),:))-0.2;
        if isempty(seq), continue; end
        for j = 1:size(seq,2)
            x = [seq(1,j),seq(2,j),seq(2,j),seq(1,j)] + IEI;
            if isWarp
                for k = 1:4
                    if x(k)>=0 && x(k)<=IEI
                       x(k) = x(k)/IEI;
                    elseif x(k)>IEI
                       x(k) = x(k)-IEI+1;
                    end
                end
            end
            y = [tickMark-0.4,tickMark-0.4,tickMark+0.4,tickMark+0.4];
            if colorByClass
                tf = vertcat(classifiedData{2:end,1})==session & vertcat(classifiedData{2:end,2})==seq(4,j);
                if contains(classifiedData{[false;tf],3},'Quality')
                    color = 'r';
                elseif contains(classifiedData{[false;tf],3},'Cue')
                    color = 'c';
                elseif contains(classifiedData{[false;tf],3},'Action')
                    color = 'b';
                else
                    color = [0.7,0.7,0.7];
                end 
            else
                color = hmmData.res.colors(seq(4,j),:);
            end
            patch(x,y,color,'edgecolor','none','facealpha',0.4);
        end
    end
    if ~isempty(trials)
        tickMark = tickMark + 2;
    end
end
handles = struct();
for i = 1:size(hmmData.res.colors,1)
    handles(i).h = plot(0.5:1,1,'color',hmmData.res.colors(i,:),'linewidth',2.5);
    handles(i).h.Visible = 'off';
end
if isWarp, xlim([-0.1,1.1]); end
xline(0,':k'); 
xline(1,':k'); 
yline(mean([max(infoMat(infoMat(:,3)==0,2)),min(infoMat(infoMat(:,3)==1,2))]),'r','linewidth',1.5);
xticks([-0.1,0,1,1.1]);
xticklabels({'($\mathrm{T}-0.1$)','(T)','(D)','($\mathrm{D}+0.1$)'});
xlabel('Warped time','interpreter','latex');
ylim([0.5,tickMark-1.5]);
myYTicks = [];
myYTickLabels = {};
for i = unique(infoMat(infoMat(:,3)==0,1)')
    block = infoMat(infoMat(:,3)==0&infoMat(:,1)==i,2);
    myYTicks = [myYTicks,mean(block)];
    myYTickLabels = [myYTickLabels,stimuli_sorted{i}(1)];
    yline(max(block)+1.5,':k');
end
yline(max(block)+11.5,':k');
for i = unique(infoMat(infoMat(:,3)==1,1)')
    block = infoMat(infoMat(:,3)==1&infoMat(:,1)==i,2);
    myYTicks = [myYTicks,mean(block)];
    myYTickLabels = [myYTickLabels,stimuli_sorted{i}(1)];
    yline(max(block)+1.5,':k');
end
yticks(myYTicks); yticklabels(myYTickLabels);
ylabel('Trials','interpreter','latex');
legend([handles(:).h],arrayfun(@(x)sprintf('%i',x),1:length(handles),'uniformoutput',false),...
    'location','eastoutside','color','none');
set(gca,'tickdir','out','color','none','box','off','ticklabelinterpreter','latex','fontsize',18);

%% Distribution of onset times of cue and action clusters
onsetTimes_cueClust = [];
onsetTimes_actClust = [];
for i = 245:254
    data = load(sprintf('%s/RawData/Simulation/simulationData%i_strongStim.mat',pwd,i));
    firingRates_thresh = data.firingRates_thresh;
    for j = 1:size(firingRates_thresh,1)
        fr_thresh_CL = firingRates_thresh{j,3}(1:60);
        fr_thresh_CR = firingRates_thresh{j,4}(1:60);
        fr_thresh_AL = firingRates_thresh{j,5}(1:60);
        fr_thresh_AR = firingRates_thresh{j,6}(1:60);
        timeCL = find(fr_thresh_CL,1)*0.05-0.025;
        timeCR = find(fr_thresh_CR,1)*0.05-0.025;
        timeAL = find(fr_thresh_AL,1)*0.05-0.025;
        timeAR = find(fr_thresh_AR,1)*0.05-0.025;
        if isempty(timeCL) && ~isempty(timeCR)
            onsetTimes_cueClust = [onsetTimes_cueClust,timeCR];
        elseif ~isempty(timeCL) && isempty(timeCR)
            onsetTimes_cueClust = [onsetTimes_cueClust,timeCL];
        elseif ~isempty(timeCL) && ~isempty(timeCR)
            onsetTimes_cueClust = [onsetTimes_cueClust,min([timeCL,timeCR])];
        end
        if isempty(timeAL) && ~isempty(timeAR)
            onsetTimes_actClust = [onsetTimes_actClust,timeAR];
        elseif ~isempty(timeAL) && isempty(timeAR)
            onsetTimes_actClust = [onsetTimes_actClust,timeAL];
        elseif ~isempty(timeAL) && ~isempty(timeAR)
            onsetTimes_actClust = [onsetTimes_actClust,min([timeAL,timeAR])];
        end
    end
end

filterConstant = 200; 
binWidth = 0.05; 
numPads = 2; 
expansion = 100; 
bins = 0:binWidth:3;

figure(1); clf; hold all;
% cue clusters
h2 = histogram(onsetTimes_cueClust,bins,'normalization','pdf','EdgeColor','none','FaceColor','k');
[~,t2,c2] = fun.hist2curve(h2,'numPads',numPads,'isPlot',false);
[t2filt,c2filt] = fun.mySmoothing(t2,c2,'zeros',0,expansion,filterConstant);
p2 = plot(t2filt,c2filt/max(c2filt),':c','linewidth',2.5);
h2.Visible = 'off';
% action clusters
h4 = histogram(onsetTimes_actClust,bins,'normalization','pdf','EdgeColor','none','FaceColor','k');
[~,t4,c4] = fun.hist2curve(h4,'numPads',numPads,'isPlot',false);
[t4filt,c4filt] = fun.mySmoothing(t4,c4,'zeros',0,expansion,filterConstant);
p4 = plot(t4filt,c4filt/max(c4filt),':b','linewidth',2.5);
h4.Visible = 'off';
% stimulus 
stim = (1/(160-150))*(exp(-(0:0.05:3000)/(160))-exp(-(0:0.05:3000)/(150)));
plot((0:0.05:3000)/1000,stim/max(stim),':k','linewidth',2.5);
% performance
for str = [25,100]
    data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr%i_silDur250_stimFall160_stimGain200.mat',pwd,str));
    perf = zeros(1,61);
    for i = 1:61
        block = vertcat(data.perf{:,i});
        perf(i) = mean(block(:,1)+0.5*(block(:,3)-block(:,2)));
    end
    plot(0:0.05:3,perf/100,'k');
end
xlim([0,3]);
xlabel('Time [s]');
ylim([0,1]);
set(gca,'tickdir','out','color','none','box','off');

%% plot any simulation trial (for Figure 3c and 4a)

% 3c (top) is Session 251, Baseline
% 3c (bottom) is Session 251, Trial 10
% 4a is Session 254, Trial 40

% -------------------------------------------------
session = 251;
trial = 10; % 'B' for baseline or # from 1-100
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
    fun.rasterPlotAbbrev(firings,data.inds,data.clustersWithRoles,totalT,net,7,true,'tick',offsetX);
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
% 4d (left) is Session 250, Trial 61 (state 11)
% 4e (left) is Session 249, Trial 72 (state 6)

% -----------------
session = 249;
trial = 72;
% -----------------

classifiedData = fun.loadVar(sprintf('%s/HMMData/Simulation/classifiedStates_sim.mat',pwd));
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
    tf = vertcat(classifiedData{2:end,1})==session & vertcat(classifiedData{2:end,2})==seq(4,i);
    if contains(classifiedData{[false;tf],3},'Quality')
        text(seq(1,i),nneurons+0.9,'Quality-coding','color','r');
    elseif contains(classifiedData{[false;tf],3},'Cue')
        text(seq(1,i),nneurons+0.9,'Cue-coding','color','c');
    elseif contains(classifiedData{[false;tf],3},'Action')
        text(seq(1,i),nneurons+0.9,'Action-coding','color','b');
    end
    patch([seq(1,i),seq(2,i),seq(2,i),seq(1,i)],...
        [0.6,0.6,nneurons+0.4,nneurons+0.4],...
        colors(seq(4,i),:),'edgecolor','none','facealpha',0.25);
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
% print(sprintf('raster_sim_hmm_S%iT%i',session,trial),'-dpdf');

%% plot any state's occurrence frequencies given trial types for any experimental session (for Figure 1b-d, right)

% 1b (right) is Session 13, State 3
% 1c (right) is Session 16, State 2
% 1d (right) is Session 2, State 2

% ----------------
session = 13;
state = 2;
% ----------------

hmmData = load(sprintf('%s/HMMData/Experiment/HMM_exp%i.mat',pwd,session));
SessionInfo = fun.loadVar(sprintf('%s/RawData/Experiment/SessionInfo.mat',pwd));
stimuli = {SessionInfo(session).TrialEvents.TrialType(:).TasteID};
%classificationData = fun.loadVar(sprintf('%s/HMMData/Experiment/classifiedStates_exp.mat',pwd));
scoredTrials = [SessionInfo(session).TrialEvents.TrialType(:).wasCorrect];
myCell = fun.getStateCell(hmmData,stimuli,scoredTrials,true);

f = figure(4); clf;
f.Position(3:4) = [332,136];
correct = [sum([myCell{[1,2],1}]==state)/sum([myCell{[1,2],2}]) sum([myCell{[3,4],1}]==state)/sum([myCell{[3,4],2}])];
incorrect = [sum([myCell{[1,2],3}]==state)/sum([myCell{[1,2],4}]) sum([myCell{[3,4],3}]==state)/sum([myCell{[3,4],4}])];
bar([correct' incorrect']);
set(gca,'XTickLabel',{'Cue Left','Cue Right'},'fontsize',10);
set(gca,'Color','none','tickdir','out','box','off');
legend('Correct Trials','Incorrect Trials','location','northeastoutside','fontsize',10);
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'renderer','painters','PaperPositionMode','auto')

f = figure(5); clf;
f.Position(3:4) = [212,136];
correct = [sum([myCell{[1,3],1}]==state)/sum([myCell{[1,3],2}]) sum([myCell{[2,4],1}]==state)/sum([myCell{[2,4],2}])];
bar(correct');
set(gca,'XTickLabel',{'Sweet','Bitter'},'fontsize',10);
set(gca,'Color','none','tickdir','out','box','off');
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'renderer','painters','PaperPositionMode','auto')

%% plot any state's occurrence frequencies given trial types for any simulation session (for Figure 4c-e, right)

% 4c (right) is Session 247, state 20
% 4d (right) is Session 250, state 11
% 4e (right) is Session 249, state 6

% ----------------
session = 249;
state = 6;
% ----------------

hmmData = load(sprintf('%s/HMMData/Simulation/HMM_sim%i.mat',pwd,session));
networkData = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,session));
%classificationData = fun.loadVar(sprintf('%s/HMMData/Simulation/classifiedStates_sim.mat',pwd));
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
set(gcf,'PaperPositionMode','auto');
% print(sprintf('stateFreq_Decision_Sess%iState%i',session,state),'-dpdf');

f = figure(7); clf;
f.Position(3:4) = [212,136];
correct = [sum([myCell{[1,3],1}]==state)/sum([myCell{[1,3],2}]) sum([myCell{[2,4],1}]==state)/sum([myCell{[2,4],2}])];
bar(correct');
set(gca,'XTickLabel',{'Sweet','Bitter'},'fontsize',10);
set(gca,'Color','none','TickLength',[0 0],'box','off');
ylabel('P(State | Trial Type)','fontsize',10);
ylim([0 1]);
set(gcf,'PaperPositionMode','auto');
% print(sprintf('stateFreq_Quality_Sess%iState%i',session,state),'-dpdf');

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
    gate.Curve = fun.getActionGate(gate.StartMean,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    
    [spkdata,~] = fun.simulation(net,neu,[],inds,[],gate.Switch,gate.Curve,S,true,warmup,BinW,timev,[],false,[]);
    
    f = figure(2); clf; hold all;
    f.Position(3:4) = [WIDTH,HEIGHT];
    firings = []; for i = 1:length(spkdata), firings = [firings; spkdata{i}]; end
    fun.rasterPlotAbbrev(firings,inds,DATA.clustersWithRoles,totalT,net,7,true,'tick',X_OFFSET);
end
