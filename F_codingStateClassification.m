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
% requires access to files:
% /RawData/Experiment/SessionInfo.mat
% /HMMData/Experiment/HMM_expX.mat for X in 1:21
% /HMMData/Experiment/HMM_expX_shuff_circ.mat for X in 1:21
% /HMMData/Experiment/HMM_expX_shuff_swap.mat for X in 1:21
% /HMMData/Simulation/HMM_simX.mat for X in 245:254
% /HMMData/Simulation/HMM_simX_shuff_circ.mat for X in 245:254
% /HMMData/Simulation/HMM_simX_shuff_swap.mat for X in 245:254
% /RawData/Simulation/simulationDataX.mat for X in 245:254
% requires access to functions:
% loadVar (in +fun)
% scoreTrials (in +fun)
% getStateCell (in +fun)
% getStateTimes (in +fun)
% classifyState (which requires chi2test and Marascuilo) (all in +fun)

%% setup
homeDir = pwd; addpath(homeDir);

%% experiment: original
classifiedStates_exp = main('experiment','',1:21,true);

%save('classifiedStates_exp.mat','classifiedStates_exp');

%% experiment: shuffled (circular)
classifiedStates_exp_shuff_circ = main('experiment','_shuff_circ',1:21,true);

%save('classifiedStates_exp_shuff_circ.mat','classifiedStates_exp_shuff_circ');

%% experiment: shuffled (swap)
classifiedStates_exp_shuff_swap = main('experiment','_shuff_swap',1:21,true);

%save('classifiedStates_exp_shuff_swap.mat','classifiedStates_exp_shuff_swap');

%% simulation: original
classifiedStates_sim = main('simulation','',245:254,false);

%save('classifiedStates_sim.mat','classifiedStates_sim');

%% simulation: shuffled (circular)
classifiedStates_sim_shuff_circ = main('simulation','_shuff_circ',245:254,false);

%save('classifiedStates_sim_shuff_circ.mat','classifiedStates_sim_shuff_circ');

%% simulation: shuffled (swap)
classifiedStates_sim_shuff_swap = main('simulation','_shuff_swap',245:254,false);

%save('classifiedStates_sim_shuff_swap.mat','classifiedStates_sim_shuff_swap');

%% main function definition
function CLASSIFIED_STATES = main(EXP_OR_SIM,MODEL,SESSIONS,ADJUST_FOR_PADDING)
CLASSIFIED_STATES = {'Session','State','classifiedLabel','onsetTimes','offsetTimes',...
                     'onsetTimes_warped','offsetTimes_warped','trialIndicator','scoredTrials'};
% setup 
counter = 1;
homeDir = pwd;
files = {}; 
if strcmp(EXP_OR_SIM,'experiment')
    cd('HMMData'); cd('Experiment'); loadDir_HMM = pwd; cd(homeDir);
    cd('RawData'); cd('Experiment'); loadDir_raw = pwd; cd(homeDir);
    for session = SESSIONS, files = [files sprintf('HMM_exp%i%s.mat',session,MODEL)]; end
elseif strcmp(EXP_OR_SIM,'simulation')
    cd('HMMData'); cd('Simulation'); loadDir_HMM = pwd; cd(homeDir);
    cd('RawData'); cd('Simulation'); loadDir_raw = pwd; cd(homeDir);
    for session = SESSIONS, files = [files sprintf('HMM_sim%i%s.mat',session,MODEL)]; end 
end
% loop over sessions
for i = 1:length(files)
    fprintf('\nAnalyzing session %i...\n',SESSIONS(i));
    % load HMM data
    cd(loadDir_HMM);
    hmmData = load(files{i});
    cd(homeDir);
    % load behavior data and score trials
    if strcmp(EXP_OR_SIM,'experiment')
        cd(loadDir_raw);
        SessionInfo = fun.loadVar('SessionInfo.mat');
        cd(homeDir);
        scoredTrials = [SessionInfo(SESSIONS(i)).TrialEvents.TrialType.wasCorrect];
        trialSequence = {SessionInfo(SESSIONS(i)).TrialEvents.TrialType.TasteID};
    elseif strcmp(EXP_OR_SIM,'simulation')
        cd(loadDir_raw);
        networkData = load(sprintf('simulationData%i.mat',SESSIONS(i)));
        cd(homeDir);
        scoredTrials = fun.scoreTrials(networkData);
        trialSequence = networkData.stimuli;
    end
    % form the stateCell for classification
    stateCell = fun.getStateCell(hmmData,trialSequence,scoredTrials,ADJUST_FOR_PADDING);
    nStates = hmmData.res.HmmParam.VarStates(hmmData.res.BestStateInd);
    trials = 1:length(scoredTrials);
    % classify each state and gather its temporal dynamics information
    for state = 1:nStates
        counter = counter + 1;
        [onsetTimes,onsetTimes_warped,...
            offsetTimes,offsetTimes_warped,trialIndicator] = fun.getStateTimes(EXP_OR_SIM,hmmData,state,trials);  
        CLASSIFIED_STATES{counter,1} = SESSIONS(i);
        CLASSIFIED_STATES{counter,2} = state;
        CLASSIFIED_STATES{counter,3} = fun.classifyState(state,stateCell,'stricterSubclassification',false);
        CLASSIFIED_STATES{counter,4} = onsetTimes;
        CLASSIFIED_STATES{counter,5} = offsetTimes;
        CLASSIFIED_STATES{counter,6} = onsetTimes_warped;
        CLASSIFIED_STATES{counter,7} = offsetTimes_warped;
        CLASSIFIED_STATES{counter,8} = trialIndicator;
        CLASSIFIED_STATES{counter,9} = scoredTrials;
    end
end
end
