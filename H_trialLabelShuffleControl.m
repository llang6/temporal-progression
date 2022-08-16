%% trialLabelShuffleControl
%
% Performs a control analysis for HMM coding state detection by comparing
% to data for which trial labels are randomly shuffled
%
% For each session, randomly permute correct and incorrect trial type
% labels. Compute state frequency by trial type, run stats to classify
% states, and see how many exclusive decision-coding and exclusive
% quality-coding states we get. Repeat and form a distribution. The number
% of coding states should be very low (chance level).
%
% This control test was only performed on the experimental data.
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Experiment/SessionInfo.mat
% /HMMData/Experiment/classifiedStates_exp.mat
% /HMMData/Experiment/HMM_expX.mat for X in 1:21
% requires access to functions:
% getStateCell (in +fun)
% classifyState (which requires chi2test and Marascuilo) (all in +fun)

%% parameters
rng(99995); % for repeatability
numIters = 10;
excludedSessions = [5,7,17,19,20];
adjustForPadding = true; % changes how states are classified by ignoring padding period if true

%% pre-load some data 
homeDir = pwd; addpath(homeDir);
fprintf('\nPre-loading data...\n');
% 'SessionInfo' (contains all behavior data)
cd('RawData'); cd('Experiment');
SessionInfo = fun.loadVar('SessionInfo.mat');
cd(homeDir);
% 'classifiedStates_exp' (original result for comparison with shuffled)
cd('HMMData'); cd('Experiment');
classifiedStates = fun.loadVar('classifiedStates_exp.mat');
cd(homeDir);
fprintf('Data pre-loaded.\n');

%% main analysis
% initialize
rng(95559); % make results repeatable
nStates_tasteID = NaN(21,numIters);
nStates_tasteQuality = NaN(21,numIters);
nStates_decision = NaN(21,numIters);
nStates_cue = NaN(21,numIters);
nStates_action = NaN(21,numIters);
nStates_non = NaN(21,numIters);
nStates_dual = NaN(21,numIters);
% loop through sessions
for session = 1:21
    if ismember(session,excludedSessions), continue; end
    fprintf('\nAnalyzing session %i...\n',session);
    % load data
    cd('HMMData'); cd('Experiment');
    hmmData = load(sprintf('HMM_exp%i.mat',session));
    cd(homeDir);
    % score trials
    scoredTrials = [SessionInfo(session).TrialEvents.TrialType.wasCorrect]==1;
    % loop through iterations
    for iter = 1:numIters
        % shuffle trial labels independently for correct and incorrect trials
        trialSequence = {SessionInfo(session).TrialEvents.TrialType.TasteID};
        trialSequence_correct = trialSequence(scoredTrials);
        trialSequence_correct_shuffled = trialSequence_correct(randperm(length(trialSequence_correct)));
        trialSequence_incorrect = trialSequence(~scoredTrials);
        trialSequence_incorrect_shuffled = trialSequence_incorrect(randperm(length(trialSequence_incorrect)));
        trialSequence_shuffled = trialSequence;
        trialSequence_shuffled(scoredTrials) = trialSequence_correct_shuffled;
        trialSequence_shuffled(~scoredTrials) = trialSequence_incorrect_shuffled;
        % create state cell for classification
        stateCell = fun.getStateCell(hmmData,trialSequence_shuffled,scoredTrials,adjustForPadding);
        % classify states and count how many of each type we found
        nStates = hmmData.res.HmmParam.VarStates(hmmData.res.BestStateInd);
        stateLabels = {};
        for state = 1:nStates
            stateLabels = [stateLabels ; fun.classifyState(state,stateCell)];
        end
        nStates_tasteID(session,iter) = sum(cellfun(@(x)strcmp(x,'Taste ID-coding'),stateLabels));
        nStates_tasteQuality(session,iter) = sum(cellfun(@(x)strcmp(x,'Exclusive Quality-coding'),stateLabels));
        nStates_decision(session,iter) = sum(cellfun(@(x)strcmp(x,'Exclusive Decision-coding'),stateLabels));
        nStates_cue(session,iter) = sum(cellfun(@(x)strcmp(x,'Cue-coding'),stateLabels));
        nStates_action(session,iter) = sum(cellfun(@(x)strcmp(x,'Action-coding'),stateLabels));
        nStates_non(session,iter) = sum(cellfun(@(x)strcmp(x,'Non-coding'),stateLabels));
        nStates_dual(session,iter) = sum(cellfun(@(x)strcmp(x,'Dual-coding'),stateLabels));
    end
end

%% plot results
x_quality_original = sum(cellfun(@(x)strcmp(x,'Exclusive Quality-coding'),classifiedStates(2:end,3))&cellfun(@(x)~ismember(x,excludedSessions),classifiedStates(2:end,1)));
x_quality_shuffled = sum(nStates_tasteQuality,'omitnan');
x_decision_original = sum(cellfun(@(x)ismember(x,{'Exclusive Decision-coding','Cue-coding','Action-coding'}),classifiedStates(2:end,3))&cellfun(@(x)~ismember(x,excludedSessions),classifiedStates(2:end,1)));
x_decision_shuffled = sum(nStates_decision,'omitnan')+sum(nStates_cue,'omitnan')+sum(nStates_action,'omitnan');
bins = min([x_quality_shuffled,x_decision_shuffled])-0.5:1:max([x_quality_shuffled,x_decision_shuffled])+0.5;
maxTick = max([x_quality_original,x_quality_shuffled,x_decision_original,x_decision_shuffled])+2;

% Quality-coding
figure(1); clf;
histogram(x_quality_shuffled,'BinEdges',bins,'FaceColor','b','EdgeColor','none'); hold on;
xline(x_quality_original,':r',{'Number found in','unshuffled data'},'LabelOrientation','horizontal');
xticks(0:maxTick);
xlim([-1,maxTick+1]);
xlabel('Number of Quality-coding states found','fontsize',14);
yticks(0:numIters);
ylim([0,numIters+1]);
ylabel('Count','fontsize',14);
title(sprintf('Effect of shuffling trial labels (%i times) on state classification',numIters),'fontsize',16);
set(gca,'TickDir','out','color','none','box','off');

% Decision-coding
figure(2); clf;
histogram(x_decision_shuffled,'BinEdges',bins,'FaceColor','b','EdgeColor','none'); hold on;
xline(x_decision_original,':r',{'Number found in','unshuffled data'},'LabelOrientation','horizontal');
xticks(0:maxTick);
xlim([-1,maxTick+1]);
xlabel('Number of Decision-coding states found','fontsize',14);
yticks(0:numIters);
ylim([0,numIters+1]);
ylabel('Count','fontsize',14);
title(sprintf('Effect of shuffling trial labels (%i times) on state classification',numIters),'fontsize',16);
set(gca,'TickDir','out','color','none','box','off');
