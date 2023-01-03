%% plotSilencing_oneNetwork
%
% Plot the data created with optoSilencing_oneNetwork.
%
% This was used to create Figure 7.
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Simulation/simulationData252_noSilencing_0_sameV0.mat
% /RawData/Simulation/simulationData252_beginningSilencing_25_sameV0.mat
% /RawData/Simulation/simulationData252_beginningSilencing_100_sameV0.mat
% /RawData/Simulation/simulationData252_cueOnsetSilencing_25_sameV0.mat
% /RawData/Simulation/simulationData252_cueOnsetSilencing_100_sameV0.mat
% /RawData/Simulation/simulationData252_middleSilencing_25_sameV0.mat
% /RawData/Simulation/simulationData252_middleSilencing_100_sameV0.mat
% /HMMData/Simulation/HMM_sim252_allConditions_multispikesremoved.mat
%
% requires access to functions:
% classifyStates (in +fun)
% chi2test (in +fun)
% Marascuilo (in +fun)
% mySmoothing (in +fun)

%% load data
all_data = struct();
data_raw_noSilencing = load(sprintf('%s/RawData/Simulation/simulationData252_noSilencing_0_sameV0.mat',pwd));
all_data(1).data = data_raw_noSilencing;
data_raw_beginningSilencing_25 = load(sprintf('%s/RawData/Simulation/simulationData252_beginningSilencing_25_sameV0.mat',pwd));
all_data(2).data = data_raw_beginningSilencing_25;
data_raw_beginningSilencing_100 = load(sprintf('%s/RawData/Simulation/simulationData252_beginningSilencing_100_sameV0.mat',pwd));
all_data(3).data = data_raw_beginningSilencing_100;
data_raw_cueSilencing_25 = load(sprintf('%s/RawData/Simulation/simulationData252_cueOnsetSilencing_25_sameV0.mat',pwd));
all_data(4).data = data_raw_cueSilencing_25;
data_raw_cueSilencing_100 = load(sprintf('%s/RawData/Simulation/simulationData252_cueOnsetSilencing_100_sameV0.mat',pwd));
all_data(5).data = data_raw_cueSilencing_100;
data_raw_middleSilencing_25 = load(sprintf('%s/RawData/Simulation/simulationData252_middleSilencing_25_sameV0.mat',pwd));
all_data(6).data = data_raw_middleSilencing_25;
data_raw_middleSilencing_100 = load(sprintf('%s/RawData/Simulation/simulationData252_middleSilencing_100_sameV0.mat',pwd));
all_data(7).data = data_raw_middleSilencing_100;
data_hmm = load(sprintf('%s/HMMData/Simulation/HMM_sim252_allConditions_multispikesremoved.mat',pwd));

choices = [];
for i = 1:length(all_data)
    choices = [choices ; all_data(i).data.choices];
end

outcomes = [];
for i = 1:length(all_data)
    outcomes0 = NaN(100,1);
    outcomes0(all_data(i).data.choices==all_data(i).data.corrChoices) = 1;
    outcomes0(all_data(i).data.choices==(~all_data(i).data.corrChoices)) = 0;
    outcomes = [outcomes ; outcomes0];
end

stimuli = repmat([repmat({'S'},25,1);repmat({'M'},25,1);repmat({'Q'},25,1);repmat({'O'},25,1)],7,1);
stimuli_fullName = repmat([repmat({'Sucrose'},25,1);repmat({'Maltose'},25,1);...
    repmat({'Quinine'},25,1);repmat({'Octaacetate'},25,1)],7,1);

classifiedStates = fun.classifyStates('simulation',252,data_hmm,outcomes,stimuli_fullName,'classificationMode',6);

%% plot all-trial decoded states across conditions
% this will also create 'stateOccurrences', which is used for 3rd plot
% --- parameters ---
colorByClass = false;
spacing_minor = 0;
spacing_major = 2;

figure(1); clf;

subplot(10,7,1);
plot((0:0.01:3),0*(0:0.01:3),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

subplot(10,7,2);
plot((0:0.01:3),25*((0:0.01:3)<=0.25),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

subplot(10,7,3);
plot((0:0.01:3),100*((0:0.01:3)<=0.25),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

subplot(10,7,4);
plot((0:0.01:3),25*((0:0.01:3)>=mean(data_raw_noSilencing.cueOnsetTimes)-0.125&(0:0.01:3)<=mean(data_raw_noSilencing.cueOnsetTimes)+0.125),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

subplot(10,7,5);
plot((0:0.01:3),100*((0:0.01:3)>=mean(data_raw_noSilencing.cueOnsetTimes)-0.125&(0:0.01:3)<=mean(data_raw_noSilencing.cueOnsetTimes)+0.125),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

subplot(10,7,6);
plot((0:0.01:3),25*((0:0.01:3)>=1.5-0.125&(0:0.01:3)<=1.5+0.125),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

subplot(10,7,7);
plot((0:0.01:3),100*((0:0.01:3)>=1.5-0.125&(0:0.01:3)<=1.5+0.125),'color',[0.9290 0.6940 0.1250],'linewidth',2);
ylim([0,100]);
set(gca,'tickdir','out','color','none','box','off');

stateOccurrences = cell(1,7);

for i = 1:7
    
    temp = false(2,100);
    stateOccurrences{i} = temp;
    
    subplot(10,7,((i+7):7:(i+63))); hold all;
    
    trials_cor = find(outcomes==1 & ismember((1:700)',((i-1)*100+1):100*i));
    trials_err = find(outcomes==0 & ismember((1:700)',((i-1)*100+1):100*i));
    trials_omi = find(isnan(outcomes) & ismember((1:700)',((i-1)*100+1):100*i));
  
    counter = 0;
    small_counter = 0;
    for b1 = 1:4
        trials_all = (((25*(b1-1)+1):25*b1) + 100*(i-1))';
        trial_blocks = {intersect(trials_all,trials_cor),intersect(trials_all,trials_err),intersect(trials_all,trials_omi)};
        for b2 = 1:length(trial_blocks)
            trial_block = trial_blocks{b2};
            stim_label_y = [];
            for j = 1:length(trial_block)
                small_counter = small_counter + 1;
                seq = data_hmm.res.hmm_postfit(trial_block(j)).sequence;
                for k = 1:size(seq,2)
                    x = [seq(1,k),seq(2,k),seq(2,k),seq(1,k)];
                    y = [j-0.4,j-0.4,j+0.4,j+0.4]+counter;
                    state_class = classifiedStates{find(vertcat(classifiedStates{2:end,2})==seq(4,k),1)+1,3};
                    if ismember(b1,[1,3]) && strcmp(state_class,'Action-coding - Left')
                        temp(1,small_counter) = true;
                    end
                    if ismember(b1,[1,3]) && strcmp(state_class,'Action-coding - Right')
                        temp(2,small_counter) = true;
                    end
                    if ismember(b1,[2,4]) && strcmp(state_class,'Action-coding - Left')
                        temp(2,small_counter) = true;
                    end
                    if ismember(b1,[2,4]) && strcmp(state_class,'Action-coding - Right')
                        temp(1,small_counter) = true;
                    end
                    if colorByClass
                        if contains(state_class,'Exclusive Quality-coding')
                            color = 'r';
                        elseif contains(state_class,'Cue-coding')
                            color = 'c';
                        elseif contains(state_class,'Action-coding')
                            color = 'b';
                        else
                            color = [0.7,0.7,0.7];
                        end
                    else
                        color = data_hmm.res.colors(seq(4,k),:);
                    end
                    patch(x,y,color,'edgecolor','none','facealpha',0.4);
                end
                stim_label_y = [stim_label_y,(y(end)+y(1))/2];
            end
            counter = counter + length(trial_block) + spacing_minor;
            if ~isempty(stim_label_y)
                yline(min(stim_label_y)-0.5,':k');
                yline(max(stim_label_y)+0.5,':k');
            end
        end
        counter = counter + spacing_major;
    end
    
    stateOccurrences{i} = temp;
    
    ylim([0,100+4*2*spacing_minor+3*spacing_major+0.6]);
    yticks([]);
    set(gca,'tickdir','out','color','none','box','off');
    
end

set(gcf,'PaperPositionMode','auto','Renderer','painters');

%print(gcf,'network252_silencingHMM_decodedStates','-depsc');

%% action coding state distribution by silencing condition
% === process data ===
myBars = NaN(3,7);
for i = 1:7
    temp = stateOccurrences{i};
    myBars(1,i) = sum(temp(1,:)&~temp(2,:));
    myBars(2,i) = sum(~temp(1,:)&temp(2,:));
    myBars(3,i) = sum(~temp(1,:)&~temp(2,:));
end

% === stats ===
cont_table = [myBars(1,:)',repmat(100,7,1)];
p = fun.chi2test(cont_table);
post_hoc = fun.Marascuilo(cont_table,0.05);

% === plot ===
figure(2); clf; hold all;
bar(myBars','stacked');
yticks(0:10:100);
xticks(1:7);
xticklabels({'Control','Weak; Beginning','Strong; Beginning','Weak; Cue onset','Strong; Cue onset',...
    'Weak; Middle','Strong; Middle'});
xlim([0,8]);
set(gca,'tickdir','out','color','none','box','off');

%% sequence comparison plots
% === process data ===
% --- parameters ---
binSize = 0.002; % [s]

% --- discretize state sequences and label all trials ---
networkSeqTrajectory = zeros(700,3/binSize);
trialLabels = cell(700,3);
t = (0:binSize:(3-binSize)) + binSize/2;

for i = 1:700
    
    seqTrajectory = 0*t;
    seq = data_hmm.res.hmm_postfit(i).sequence;
    for j = 1:size(seq,2)
        seqTrajectory(t>=seq(1,j)&t<seq(2,j)) = seq(4,j);
    end
    networkSeqTrajectory(i,:) = seqTrajectory;
    
    network_condition = ceil(i/100);
    switch network_condition
        case 1
            trialLabels{i,1} = 'control';
        case 2
            trialLabels{i,1} = 'beg25';
        case 3
            trialLabels{i,1} = 'beg100';
        case 4
            trialLabels{i,1} = 'cue25';
        case 5
            trialLabels{i,1} = 'cue100';
        case 6
            trialLabels{i,1} = 'mid25';
        case 7
            trialLabels{i,1} = 'mid100';
    end
    
    trial = i - 100*(network_condition - 1);
    switch ceil(trial/25)
        case 1
            trialLabels{i,2} = 'S';
        case 2
            trialLabels{i,2} = 'M';
        case 3
            trialLabels{i,2} = 'Q';
        case 4
            trialLabels{i,2} = 'O';
    end
    
    choice = all_data(network_condition).data.choices(trial);
    corrChoice = all_data(network_condition).data.corrChoices(trial);
    if choice == corrChoice
        trialLabels{i,3} = 'correct';
    elseif choice == ~corrChoice
        trialLabels{i,3} = 'error';
    else
        trialLabels{i,3} = 'omitted';
    end
    
end

% --- compute state sequence similarities ---
% control vs. control
tr1 = getTrials(trialLabels,'control','S','all'); 
tr2 = getTrials(trialLabels,'control','S','all');
similarity.controlVsControl_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all'); 
tr2 = getTrials(trialLabels,'control','M','all');
similarity.controlVsControl_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all'); 
tr2 = getTrials(trialLabels,'control','Q','all');
similarity.controlVsControl_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all'); 
tr2 = getTrials(trialLabels,'control','O','all');
similarity.controlVsControl_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsControl = mean([similarity.controlVsControl_S;...
    similarity.controlVsControl_M;...
    similarity.controlVsControl_Q;...
    similarity.controlVsControl_O],1);

% control vs. beginning 25
tr1 = getTrials(trialLabels,'control','S','all');
tr2 = getTrials(trialLabels,'beg25','S','all');
similarity.controlVsBeg25_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all');
tr2 = getTrials(trialLabels,'beg25','M','all');
similarity.controlVsBeg25_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all');
tr2 = getTrials(trialLabels,'beg25','Q','all');
similarity.controlVsBeg25_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all');
tr2 = getTrials(trialLabels,'beg25','O','all');
similarity.controlVsBeg25_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsBeg25 = mean([similarity.controlVsBeg25_S;...
    similarity.controlVsBeg25_M;...
    similarity.controlVsBeg25_Q;...
    similarity.controlVsBeg25_O],1);

% control vs. beginning 100
tr1 = getTrials(trialLabels,'control','S','all');
tr2 = getTrials(trialLabels,'beg100','S','all');
similarity.controlVsBeg100_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all');
tr2 = getTrials(trialLabels,'beg100','M','all');
similarity.controlVsBeg100_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all');
tr2 = getTrials(trialLabels,'beg100','Q','all');
similarity.controlVsBeg100_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all');
tr2 = getTrials(trialLabels,'beg100','O','all');
similarity.controlVsBeg100_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsBeg100 = mean([similarity.controlVsBeg100_S;...
    similarity.controlVsBeg100_M;...
    similarity.controlVsBeg100_Q;...
    similarity.controlVsBeg100_O],1);

% control vs. cue 25
tr1 = getTrials(trialLabels,'control','S','all');
tr2 = getTrials(trialLabels,'cue25','S','all');
similarity.controlVsCue25_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all');
tr2 = getTrials(trialLabels,'cue25','M','all');
similarity.controlVsCue25_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all');
tr2 = getTrials(trialLabels,'cue25','Q','all');
similarity.controlVsCue25_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all');
tr2 = getTrials(trialLabels,'cue25','O','all');
similarity.controlVsCue25_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsCue25 = mean([similarity.controlVsCue25_S;...
    similarity.controlVsCue25_M;...
    similarity.controlVsCue25_Q;...
    similarity.controlVsCue25_O],1);

% control vs. cue 100
tr1 = getTrials(trialLabels,'control','S','all');
tr2 = getTrials(trialLabels,'cue100','S','all');
similarity.controlVsCue100_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all');
tr2 = getTrials(trialLabels,'cue100','M','all');
similarity.controlVsCue100_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all');
tr2 = getTrials(trialLabels,'cue100','Q','all');
similarity.controlVsCue100_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all');
tr2 = getTrials(trialLabels,'cue100','O','all');
similarity.controlVsCue100_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsCue100 = mean([similarity.controlVsCue100_S;...
    similarity.controlVsCue100_M;...
    similarity.controlVsCue100_Q;...
    similarity.controlVsCue100_O],1);

% control vs. mid 25
tr1 = getTrials(trialLabels,'control','S','all');
tr2 = getTrials(trialLabels,'mid25','S','all');
similarity.controlVsMid25_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all');
tr2 = getTrials(trialLabels,'mid25','M','all');
similarity.controlVsMid25_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all');
tr2 = getTrials(trialLabels,'mid25','Q','all');
similarity.controlVsMid25_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all');
tr2 = getTrials(trialLabels,'mid25','O','all');
similarity.controlVsMid25_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsMid25 = mean([similarity.controlVsMid25_S;...
    similarity.controlVsMid25_M;...
    similarity.controlVsMid25_Q;...
    similarity.controlVsMid25_O],1);

% control vs. mid 100
tr1 = getTrials(trialLabels,'control','S','all');
tr2 = getTrials(trialLabels,'mid100','S','all');
similarity.controlVsMid100_S = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','M','all');
tr2 = getTrials(trialLabels,'mid100','M','all');
similarity.controlVsMid100_M = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','Q','all');
tr2 = getTrials(trialLabels,'mid100','Q','all');
similarity.controlVsMid100_Q = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
tr1 = getTrials(trialLabels,'control','O','all');
tr2 = getTrials(trialLabels,'mid100','O','all');
similarity.controlVsMid100_O = sim_t(networkSeqTrajectory(tr1,:),networkSeqTrajectory(tr2,:));
similarity.controlVsMid100 = mean([similarity.controlVsMid100_S;...
    similarity.controlVsMid100_M;...
    similarity.controlVsMid100_Q;...
    similarity.controlVsMid100_O],1);

% === plot ===
% --- parameters ---
filter_constant = 150; % [bins]
fontsize_small = 12;
fontsize_big = 12;

figure(3); clf; hold all;
yline(1,':k');
y = similarity.controlVsBeg25./similarity.controlVsControl;
[t_filt,y_filt] = fun.mySmoothing(t,y,'continuation',100,1,filter_constant);
p1 = plot(t_filt,y_filt,':g','linewidth',2);
y = similarity.controlVsBeg100./similarity.controlVsControl;
[t_filt,y_filt] = fun.mySmoothing(t,y,'continuation',100,1,filter_constant);
p2 = plot(t_filt,y_filt,'g','linewidth',2);
y = similarity.controlVsCue25./similarity.controlVsControl;
[t_filt,y_filt] = fun.mySmoothing(t,y,'continuation',100,1,filter_constant);
p3 = plot(t_filt,y_filt,':b','linewidth',2);
y = similarity.controlVsCue100./similarity.controlVsControl;
[t_filt,y_filt] = fun.mySmoothing(t,y,'continuation',100,1,filter_constant);
p4 = plot(t_filt,y_filt,'b','linewidth',2);
y = similarity.controlVsMid25./similarity.controlVsControl;
[t_filt,y_filt] = fun.mySmoothing(t,y,'continuation',100,1,filter_constant);
p5 = plot(t_filt,y_filt,':r','linewidth',2);
y = similarity.controlVsMid100./similarity.controlVsControl;
[t_filt,y_filt] = fun.mySmoothing(t,y,'continuation',100,1,filter_constant);
p6 = plot(t_filt,y_filt,'r','linewidth',2);
xlim([0,3]);
legend([p1,p2,p3,p4,p5,p6],{'Beginning (25%)',...
    'Beginning (100%)',...
    'Cue (25%)',...
    'Cue (100%)',...
    'Middle (25%)',...
    'Middle (100%)'},'location','northeastoutside');
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);
ylabel('Sequence similarity relative to control','fontsize',fontsize_big);

set(gcf,'renderer','painters','paperpositionmode','auto');

%print(gcf,'sequenceSimilarities','-dpdf');

%% local function definitions
function SIMILARITY = sim_t(MAT1,MAT2)
% the measure of similarity between two vectors of states over trials used
% here is equivalent to the inner product of the estimated probability
% distributions over states
% sim(x,x0) = <p_x,p_x0>
if size(MAT1,2)~=size(MAT2,2)
    error('Input matrices must have the same size in the second dimension');
end
SIMILARITY = NaN(1,size(MAT1,2));
for t = 1:size(MAT1,2)
    temp = NaN(size(MAT1,1),size(MAT2,1));
    for i = 1:size(MAT1,1)
        for j = 1:size(MAT2,1)
            temp(i,j) = MAT1(i,t)==MAT2(j,t);
        end
    end
    SIMILARITY(t) = mean(temp,'all');
end
end

function TRIALS = getTrials(LABEL_CELL,C1,C2,C3)
if strcmp(C1,'all')
    c1 = true(size(LABEL_CELL,1),1);
else
    c1 = cellfun(@(x)strcmp(x,C1),LABEL_CELL(:,1));
end
if strcmp(C2,'all')
    c2 = true(size(LABEL_CELL,1),1);
else
    c2 = cellfun(@(x)strcmp(x,C2),LABEL_CELL(:,2));
end
if strcmp(C3,'all')
    c3 = true(size(LABEL_CELL,1),1);
else
    c3 = cellfun(@(x)strcmp(x,C3),LABEL_CELL(:,3));
end
TRIALS = find(c1&c2&c3);
end
