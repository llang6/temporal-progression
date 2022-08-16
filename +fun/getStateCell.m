function STATE_CELL = getStateCell(HMMDATA,TRIAL_SEQUENCE,SCORED_TRIALS,ADJUST_FOR_PADDING)
%STATE_CELL = getStateCell(HMMDATA,TRIAL_SEQUENCE,SCORED_TRIALS,ADJUST_FOR_PADDING)
%
% Creates the state cell object (a contingency table of sorts) for state
% classification
%
% Inputs:
% HMMDATA: a structure of loaded HMM data, e.g. HMMDATA = load('HMM_exp1.mat');
% TRIAL_SEQUENCE: a cell containing the identities of presented stimuli (as
% strings) on each trial
% SCORED_TRIALS: a column vector of 0s and 1s, 0 for incorrect trial and 1
% for correct trials
% ADJUST_FOR_PADDING: true/false flag; when dealing with experimental data,
% HMMs were fit to a window 0.1 s before and after the taste/decision
% interevent interval; if this flag is set true, state occurrences
% *entirely contained* within these 0.1 s edge padding intervals will *not*
% be considered occurrences of the states
%
% Output:
% STATE_CELL: a 4x4 cell; each row corresponds to a stimulus identity
% ('Sucrose','Quinine','Maltose','Octaacetate' in that order); column 1
% contains vectors of IDs of states that occurred across all trials of the
% given type; column 2 contains a single number, the total number of trials
% of the given type; columns 3 and 4 mirror columns 1 and 2 but columns 1
% and 2 are for correct trials whereas columns 3 and 4 are for incorrect
% trials
%
% Occurrence frequencies of states in any given trial type are easily
% calculated from STATE_CELL, e.g. occurrence frequency of state X in
% trials of type...
% Correct   'Sucrose': sum(STATE_CELL{1,1}==X)/STATE_CELL{1,2}
% Incorrect 'Sucrose': sum(STATE_CELL{1,3}==X)/STATE_CELL{1,4}
% Correct   'Maltose': sum(STATE_CELL{3,1}==X)/STATE_CELL{3,2}
% Incorrect 'Maltose': sum(STATE_CELL{3,3}==X)/STATE_CELL{3,4}
%
% -LL
%
tastes = {'Sucrose','Quinine','Maltose','Octaacetate'};
STATE_CELL = cell(4,4);
for j = 1:4
    trialSeq_correct = find(reshape(strcmp(TRIAL_SEQUENCE,tastes{j}),1,[])&reshape(SCORED_TRIALS==1,1,[]));
    trialSeq_incorrect = find(reshape(strcmp(TRIAL_SEQUENCE,tastes{j}),1,[])&reshape(SCORED_TRIALS==0,1,[]));
    states_allTrials_correct = [];
    states_allTrials_incorrect = [];
    for k = 1:length(trialSeq_correct)
        % count state occurrences in correct trials
        states_trial = [];
        state_info = HMMDATA.res.hmm_postfit(trialSeq_correct(k)).sequence;
        if ~isempty(state_info)
            if ADJUST_FOR_PADDING
                % accept no states with offsets before left padding ends
                ind1 = state_info(2,:) > (HMMDATA.win_train(trialSeq_correct(k),1) + 0.1);
                % or onsets after right padding begins
                ind2 = state_info(1,:) < 0;
                states_trial = [states_trial, state_info(4,(ind1&ind2))];
            else
                states_trial = [states_trial, state_info(4,:)];
            end
        end
        states_allTrials_correct = [states_allTrials_correct, unique(states_trial)];
    end
    for k = 1:length(trialSeq_incorrect)
        % count state occurrences in incorrect trials
        states_trial = [];
        state_info = HMMDATA.res.hmm_postfit(trialSeq_incorrect(k)).sequence;
        if ~isempty(state_info)
            if ADJUST_FOR_PADDING
                % accept no states with offsets before left padding ends
                ind1 = state_info(2,:) > (HMMDATA.win_train(trialSeq_incorrect(k),1) + 0.1);
                % or onsets after right padding begins
                ind2 = state_info(1,:) < 0;
                states_trial = [states_trial, state_info(4,(ind1&ind2))];
            else
                states_trial = [states_trial, state_info(4,:)];
            end
        end
        states_allTrials_incorrect = [states_allTrials_incorrect, unique(states_trial)];
    end
    % Error check %%%%%%%%%%%%%%%%%%%%
    for n = unique(states_allTrials_correct); if sum(states_allTrials_correct==n) > length(trialSeq_correct); error('Check the state-finding algorithm -- there should not be more of any state than total trials'); end; end
    for n = unique(states_allTrials_incorrect); if sum(states_allTrials_incorrect==n) > length(trialSeq_incorrect); error('Check the state-finding algorithm -- there should not be more of any state than total trials'); end; end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STATE_CELL{j,1} = states_allTrials_correct;
    STATE_CELL{j,2} = length(trialSeq_correct);
    STATE_CELL{j,3} = states_allTrials_incorrect;
    STATE_CELL{j,4} = length(trialSeq_incorrect);
end
end