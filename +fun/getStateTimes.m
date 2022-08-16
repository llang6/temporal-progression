function [ONSET_TIMES,ONSET_TIMES_WARPED,...
    OFFSET_TIMES,OFFSET_TIMES_WARPED,TRIAL_INDICATOR] = getStateTimes(EXP_OR_SIM,HMMDATA,STATE,TRIALS)
%[ONSET_TIMES,ONSET_TIMES_WARPED,...
% OFFSET_TIMES,OFFSET_TIMES_WARPED,TRIAL_INDICATOR] = getStateTimes(EXP_OR_SIM,HMMDATA,STATE,TRIALS)
%
% Extract all within-trial state occurrence times for a given state
%
% Inputs:
% EXP_OR_SIM: 'experiment' or 'simulation' as a string; indicates whether HMM data comes
% from experimental or simulated data
% HMMDATA: a structure of loaded HMM data, e.g. HMMDATA = load('HMM_exp1.mat');
% STATE: the ID (a single number) of the state of interest
% TRIALS: a row vector of all the trials to consider when extracting
% occurrence times
%
% Outputs:
% ONSET_TIMES: vector of all within-trial state onset times
% ONSET_TIMES_WARPED: above times divided by the taste/decision interevent
% interval
% OFFSET_TIMES: vector of all within-trial state offset times
% OFFSET_TIMES_WARPED: above times divided by the taste/decision interevent
% interval 
% TRIAL_INDICATOR: vector of trials from which each onset/offset time was
% extracted (ONSET_TIMES, OFFSET_TIMES, and TRIAL_INDICATOR all have the
% same length)
%
% -LL
%
ONSET_TIMES = []; 
ONSET_TIMES_WARPED = []; 
OFFSET_TIMES = []; 
OFFSET_TIMES_WARPED = [];
TRIAL_INDICATOR = [];
for trial = TRIALS
    sequence = HMMDATA.res.hmm_postfit(trial).sequence; 
    if isempty(sequence), continue; end 
    statesOccurred = sequence(4,:); 
    if strcmpi(EXP_OR_SIM,'experiment')
        % correcting for padding in experiment HMM fitting
        warp = diff(HMMDATA.win_train(trial,:)) - 0.2; % equivalent to but more generalizable than below
        %warp = abs(HMMDATA.win_train(trial,1) + 0.1);
    elseif strcmpi(EXP_OR_SIM,'simulation')
        warp = diff(HMMDATA.win_train(trial,:));
        %warp = 3;
    else
        error('Invalid input for ''EXP_OR_SIM''. Valid inputs are ''experiment'' or ''simulation''.');
    end
    %newOnsets = sequence(1,find(statesOccurred==state,1)); % only count first onset time
    newOnsets = sequence(1,statesOccurred==STATE); % count all onset times
    %newOnsets = newOnsets - (HMMDATA.win_train(trial,1) + 0.1); % relative to taste
    newOnsets_warped = newOnsets/warp; % time "warping"
    newOffsets = sequence(2,statesOccurred==STATE);
    newOffsets_warped = newOffsets/warp;
    ONSET_TIMES = [ONSET_TIMES, newOnsets]; 
    ONSET_TIMES_WARPED = [ONSET_TIMES_WARPED, newOnsets_warped];
    OFFSET_TIMES = [OFFSET_TIMES, newOffsets]; 
    OFFSET_TIMES_WARPED = [OFFSET_TIMES_WARPED, newOffsets_warped];
    TRIAL_INDICATOR = [TRIAL_INDICATOR, repmat(trial,1,length(newOnsets))];
end  
end