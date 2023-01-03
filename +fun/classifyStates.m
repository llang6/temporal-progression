function CLASSIFIED_STATES = classifyStates(EXP_OR_SIM,SESSION,HMM_DATA,SCORED_TRIALS,STIMULI,varargin)
%CLASSIFIED_STATES = classifyStates(EXP_OR_SIM,SESSION,HMM_DATA,SCORED_TRIALS,STIMULI) 
%
%Classifies all the states contained in HMM_DATA according to differential
%occurrence frequencies in trial types defined by SCORED_TRIALS and STIMULI
%
% This function is rather specific and not very flexible. Follow this
% format:
%
% Input EXP_OR_SIM: a string, either 'experiment' or 'simulation' (lets the
% function know whether it should adjust for padding, which only applies to
% HMMs fit to experimental data)
%
% Input SESSION: an integer denoting the session number (simply used as a
% label in the output)
%
% Input HMM_DATA: a structure containing all the HMM data
%                 ex.
%                 hmm_data = load(sprintf('%s/HMMData/Experiment/HMM_exp1.mat',pwd));
%                 classified_states = classifyStates('experiment',1,hmm_data,...);
%
% Input SCORED_TRIALS: a column vector with 1 meaning correct trial, 0
% meaning incorrect trial, and NaN meaning omitted trial
%
% Input STIMULI: a cell array with values 'Sucrose', 'Quinine', 'Maltose',
% or 'Octaacetate' to indicate the stimulus given on each trial
%
% Output CLASSIFIED_STATES: a cell with all state classification
% information (see the headers in the output)
%
% Optional 'Name',Value pair arguments can be provided after STIMULI:
%     'isHeader': logical (default true), controls whether the output CLASSIFIED_STATES
%     cell will have a labeled header row
%
%     'classificationMode': an integer from 1 to 6 (default 3), controls the
%     classification mode used by the function classifyState (see
%     classifyState for information regarding each mode)
%
%     'stricterSubclass': logical (default false), passed to classifyState
%     to control whether an additional classification criterion
%     (Decision-coding states need at least 5 incorrectly performed trials
%     to be sub-classified into Cue- or Action-coding) should be applied
%
% -LL
%

% default parameters
default_isHeader = true;
default_classificationMode = 3;
default_stricterSubclass = false;
% parse input
if isempty(varargin)
    isHeader = default_isHeader;
    classificationMode = default_classificationMode;
    stricterSubclass = default_stricterSubclass;
else
    ind = find(cellfun(@(x)strcmpi('isHeader',x),varargin));
    if isempty(ind), isHeader = default_isHeader; else, isHeader = varargin{ind+1}; end
    ind = find(cellfun(@(x)strcmpi('classificationMode',x),varargin));
    if isempty(ind), classificationMode = default_classificationMode; else, classificationMode = varargin{ind+1}; end
    ind = find(cellfun(@(x)strcmpi('stricterSubclassification',x),varargin));
    if isempty(ind), stricterSubclass = default_stricterSubclass; else, stricterSubclass = varargin{ind+1}; end
end
if strcmp(EXP_OR_SIM,'experiment')
    adjust_for_padding = true;
elseif strcmp(EXP_OR_SIM,'simulation')
    adjust_for_padding = false;
end
% initialize output data
if isHeader
    CLASSIFIED_STATES = {'Session','State','classifiedLabel','onsetTimes','offsetTimes',...
        'onsetTimes_warped','offsetTimes_warped','trialIndicator','scoredTrials'};
    counter = 1;
else
    CLASSIFIED_STATES = {};
    counter = 0;
end
% form the stateCell for classification
stateCell = fun.getStateCell(HMM_DATA,STIMULI,SCORED_TRIALS,adjust_for_padding);
nStates = HMM_DATA.res.HmmParam.VarStates(HMM_DATA.res.BestStateInd);
trials = 1:length(SCORED_TRIALS);
% classify each state and gather its temporal dynamics information
for state = 1:nStates
    counter = counter + 1;
    [onsetTimes,onsetTimes_warped,...
        offsetTimes,offsetTimes_warped,trialIndicator] = fun.getStateTimes(EXP_OR_SIM,HMM_DATA,state,trials);  
    CLASSIFIED_STATES{counter,1} = SESSION;
    CLASSIFIED_STATES{counter,2} = state;
    CLASSIFIED_STATES{counter,3} = fun.classifyState(state,stateCell,'classificationMode',classificationMode,...
        'stricterSubclassification',stricterSubclass);
    CLASSIFIED_STATES{counter,4} = onsetTimes;
    CLASSIFIED_STATES{counter,5} = offsetTimes;
    CLASSIFIED_STATES{counter,6} = onsetTimes_warped;
    CLASSIFIED_STATES{counter,7} = offsetTimes_warped;
    CLASSIFIED_STATES{counter,8} = trialIndicator;
    CLASSIFIED_STATES{counter,9} = SCORED_TRIALS;
end
end