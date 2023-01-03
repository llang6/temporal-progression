%% codingStateClassification
% 
% Classify all decoded HMM states and package each labeled state along with
% its within-trial onset times into an organized cell called
% 'classifiedStates'.
%
% This script generates the 'classifiedStates' cell for all data types
% (original, circular shuffle, swap shuffle) and for both experiment and
% simulation data.
%
% You can also just execute code blocks individually to generate a
% particular 'classifiedStates' cell.
%
% Essentially the script loads HMM data for state information (which states
% appeared in which trials and when) and raw data for trial information
% (what was the stimulus for each trial, was it performed correctly, and
% when were the taste/decision events).
% It organizes all of this into an intermediate piece of data called
% 'stateCell,' then passes this to the classifyState function to perform
% the statistical analysis, which involves Chi-squared tests run on
% contingency tables derived from 'stateCell' (see classifyState in +fun
% for more information).
% In addition, each state's within-trial onset and offset times are
% collected. The times are saved in their original format (relative to
% either the taste/decision event) and in a warped format where a unit of
% time is the taste/decision interevent interval. The warped format allows
% for comparison across trials with different interevent intervals and for
% easy conversion between times relative-to-taste and relative-to-decision.
%
% The final output cell 'classifiedStates' is now all you need to plot
% information about the temporal dynamics of the coding states.
%
% - LL
%

%% dependencies
% requires access to data:
% /RawData/Experiment/SessionInfo.mat
% /HMMData/Experiment/HMM_expX.mat for X in 1:21
% /HMMData/Experiment/HMM_expX_shuff_circ.mat for X in 1:21
% /HMMData/Experiment/HMM_expX_shuff_swap.mat for X in 1:21
% /HMMData/Simulation/HMM_simX.mat for X in 245:254
% /HMMData/Simulation/HMM_simX_shuff_circ.mat for X in 245:254
% /HMMData/Simulation/HMM_simX_shuff_swap.mat for X in 245:254
% /RawData/Simulation/simulationDataX.mat for X in 245:254
%
% requires access to functions:
% loadVar (in +fun)
% scoreTrials (in +fun)
% classifyStates (which requires other functions, all in +fun)

%% setup
addpath(pwd);

%% experiment: unshuffled
classifiedStates_exp = main('experiment','','',1:21,false);
%save(sprintf('%s/HMMData/Experiment/classifiedStates_exp.mat',pwd),'classifiedStates_exp');

%% experiment: unshuffled (stricter sub-classification)
classifiedStates_exp_stricterSubClassification = main('experiment','','',1:21,true);
%save(sprintf('%s/HMMData/Experiment/classifiedStates_exp_stricterSubClassification.mat',pwd),'classifiedStates_exp_stricterSubClassification');

%% experiment: shuffled (circular)
classifiedStates_exp_shuff_circ = main('experiment','_shuff_circ','',1:21,false);
%save(sprintf('%s/HMMData/Experiment/classifiedStates_exp_shuff_circ.mat',pwd),'classifiedStates_exp_shuff_circ');

%% experiment: shuffled (swap)
classifiedStates_exp_shuff_swap = main('experiment','_shuff_swap','',1:21,false);
%save(sprintf('%s/HMMData/Experiment/classifiedStates_exp_shuff_swap.mat',pwd),'classifiedStates_exp_shuff_swap');

%% simulation: unshuffled, original stim
classifiedStates_sim = main('simulation','','',245:254,false);
%save(sprintf('%s/HMMData/Simulation/classifiedStates_sim.mat',pwd),'classifiedStates_sim');

%% simulation: unshuffled, long stim
classifiedStates_sim = main('simulation','','_longStim',245:254,false);
%save(sprintf('%s/HMMData/Simulation/classifiedStates_sim_longStim.mat',pwd),'classifiedStates_sim');

%% simulation: unshuffled, strong stim
classifiedStates_sim = main('simulation','','_strongStim',245:254,false);
%save(sprintf('%s/HMMData/Simulation/classifiedStates_sim_strongStim.mat',pwd),'classifiedStates_sim');

%% simulation: shuffled (circular)
classifiedStates_sim_shuff_circ = main('simulation','_shuff_circ','',245:254,false);
%save(sprintf('%s/HMMData/Simulation/classifiedStates_sim_shuff_circ.mat',pwd),'classifiedStates_sim_shuff_circ');

%% simulation: shuffled (swap)
classifiedStates_sim_shuff_swap = main('simulation','_shuff_swap','',245:254,false);
%save(sprintf('%s/HMMData/Simulation/classifiedStates_sim_shuff_swap.mat',pwd),'classifiedStates_sim_shuff_swap');

%% main function definition
function CLASSIFIED_STATES = main(EXP_OR_SIM,MODEL,STIM_TYPE,SESSIONS,STRICTER_SUBCLASS)
% setup 
files = {}; 
if strcmp(EXP_OR_SIM,'experiment')
    for session = SESSIONS
        files = [files, sprintf('%s/HMMData/Experiment/HMM_exp%i%s.mat',pwd,session,MODEL)]; 
    end
elseif strcmp(EXP_OR_SIM,'simulation')
    for session = SESSIONS 
        files = [files, sprintf('%s/HMMData/Simulation/HMM_sim%i%s%s.mat',pwd,session,MODEL,STIM_TYPE)]; 
    end 
end
% loop over sessions
CLASSIFIED_STATES = {};
for i = 1:length(files)
    fprintf('\nAnalyzing session %i...\n',SESSIONS(i));
    % load HMM data
    hmmData = load(files{i});
    % load behavior data and score trials
    if strcmp(EXP_OR_SIM,'experiment')
        SessionInfo = fun.loadVar(sprintf('%s/RawData/Experiment/SessionInfo.mat',pwd));
        scoredTrials = [SessionInfo(SESSIONS(i)).TrialEvents.TrialType.wasCorrect];
        trialSequence = {SessionInfo(SESSIONS(i)).TrialEvents.TrialType.TasteID};
    elseif strcmp(EXP_OR_SIM,'simulation')
        networkData = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,SESSIONS(i)));
        scoredTrials = fun.scoreTrials(networkData);
        trialSequence = networkData.stimuli;
    end
    if i == 1, isHeader = true; else, isHeader = false; end
    CLASSIFIED_STATES = [CLASSIFIED_STATES ; fun.classifyStates(EXP_OR_SIM,SESSIONS(i),hmmData,...
        scoredTrials,trialSequence,'stricterSubclass',STRICTER_SUBCLASS,'isHeader',isHeader)];
end
end
