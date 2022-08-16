%% optoSilencingExperiment
%
% This script carries out the simulated optogenetic perturbation
% experiments.
% 
% Silencing is a square pulse input to all inhibitory neurons in the
% network. Silencing is parameterized by its strength ("silence_strength,"
% the % increase in baseline external current), its duration
% ("silence_duration"; 250 ms used for all simulations in the paper), and
% its center value ("silence_center").
%
% %%%%%%%%
% PART 1
% %%%%%%%%
%
% Effect of silencing and stimulus parameters on network task accuracy
% 
% The code here allows one to define the parameters of the silencing. 
% It also allows one to parameterize the taste stimulus presented to the
% network by changing its decay time constant ("stim_tauFall"; 160 ms is
% the baseline/standard).
%
% For any given parameter set {silence_strength, silence_duration,
% stim_tauFall}, this script will calculate network task accuracy as a
% function of silence_center by sliding the center value of the silencing
% across the 3 s trial window and, at each point, simulating 1000 trials
% (10 networks, 100 trials each).
% It, therefore, generates the data plotted in a single curve in Figure 6C
% or Figure 6D (middle).
% Results are saved in a session/network-by-center cell called "perf."
% Each entry in the cell summarizes the performance of each network for
% parameterized silencing centered at each point.
% This summary is in the form of a 3D vector with entries: number of
% correct trials, number of non-omitted trials, number of total trials.
%
% This is simulation-intensive, since 1000 trials are simulated per center
% value (and there are 61 center values if you slide over the entire trial
% window). 
% The "for loop" over sessions (networks) was changed to a "parfor loop" to
% generate the data for the paper.
%
% %%%%%%%%
% PART 2
% %%%%%%%%
%
% Effect of silencing strength on cluster firing rates
%
% This part of the script generates the data shown in Figure 6B. 
% For a given {silence_strength, silence_duration, silence_center}, the
% code will simply simulate 1000 trials (again, 10 networks with 100 trials
% each) and compute pre-silencing and during-silencing firing rates for
% each cluster in each trial.
% Results are saved in a 10-dimensional structure (1D per session/network),
% with each dimension having 4 different fields, which are cluster-by-trial
% firing rate matrices (1 per combination of excitatory/inhibitory clusters
% before/during silencing).
% Theoretically the data generated by these simulations could have been
% acquired from the simulations run in PART 1, but in practice they were
% generated separately.
% The "for loop" over sessions/networks was also changed to a "parfor loop"
% before running.
% 
% %%%%%%%%
% PART 3
% %%%%%%%%
%
% Plot a single trial simulation
%
% This section generates Figure 6A. 
% For a given network/session, trial, silencing parameterization, and
% stimulus parameterization, this code simulates the trial with silencing
% and plots the external inputs, spike raster, and cluster firing rate
% curves.
% 
% %%%%%%%%
% PART 4
% %%%%%%%%
%
% Generate Figure plots 
%
% This code loads the data generated and saved by Parts 1-2 (as well as
% some additional baseline data) and creates the plots from Figure 6B-D.
% 
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Simulation/simulationDataX.mat for X in 245:254
% /HMMData/Simulation/onsetTimeDist_XCS.mat for X in {Q,A,C}
% /ProcessedData/Simulation/deltaFRs_strX_durY_cZ.mat for various X, Y, Z
% /ProcessedData/Simulation/perf_silStrX_silDurY_stimFallZ_stimGain60.mat for various X, Y, Z
% requires access to functions:
% gaussfilt (in +fun)
% distinguishable_colors (in +fun)

%% PART 1

% options and global parameters -------------------------
sessions = 1:10;
numTrials = 100;
warmup = 0; % simulate for this long before considering the trial to have begun, [ms]
totalT = 3000; % trial time, [ms]
dt = 0.05; % [ms]

silence_center = 0:0.05:3; % relative to stimulus onset, [s]
silence_duration = 0.25; % [s]
silence_strength = 100; % [% increase]

changeStim = true;
stim_Gain = 0.60; % [fraction increase]
stim_tauRise = 0.150; % [s]
stim_tauFall = 0.705; % [s]

autoSaveOutput = false;

% main simulations -------------------------
timev = 0:dt:warmup+totalT;
index2 = 1:length(silence_center);
stim_Curve(timev>=warmup) = (1/((stim_tauFall*1000)-(stim_tauRise*1000)))*...
    (exp(-(0:dt:totalT)/(stim_tauFall*1000))-...
    exp(-(0:dt:totalT)/(stim_tauRise*1000)));
stim_Curve = stim_Curve/max(abs(stim_Curve))*stim_Gain;
perf = cell(10,length(silence_center));

for session_count = sessions

    % session-specific parameters
    data = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,244+session_count));
    BinW = data.BinW;
    binSize = data.binSize;
    gate_startTimes = data.gate.StartTimes;
    gate_duration = data.gate.Duration;
    gate_floor = data.gate.Floor;
    gate_ceiling = data.gate.Ceiling;
    gate_switch = data.gate.Switch;
    net = data.net;
    neu = data.neu;
    stim = data.stim;
    if changeStim, stim.Curve = stim_Curve; end
    inds = data.inds;
    S = data.S;
    stateThresh = data.stateThresh;
    stateMinTime = data.stateMinTime;
    bound1 = data.bound1;
    bound2 = data.bound2;

    for center_count = index2

        silence_curve = getSilencing(silence_center(center_count),silence_duration,silence_strength,warmup,totalT,dt);

        fprintf('\n\nSession %i, silencing centered at %.3f:\n\n',244+session_count,silence_center(center_count));

        firings_all = cell(numTrials,1);
        firingRates = cell(numTrials,6);
        stimuli = [repmat({'Sucrose'},1,ceil(numTrials/4)),... 
                   repmat({'Maltose'},1,ceil(numTrials/4)),...
                   repmat({'Quinine'},1,ceil(numTrials/4)),...
                   repmat({'Octaacetate'},1,ceil(numTrials/4))];

        for iter = 1:numTrials

            timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\nSimulating trial %i...',iter);

            gate_curve = getActionGate(gate_startTimes(iter),gate_duration,gate_floor,gate_ceiling,warmup,totalT,dt);
            [spkdata,~] = simulation(net,neu,stim,inds,stimuli{iter},gate_switch,gate_curve,S,true,...
                warmup,BinW,timev,silence_curve,false);

            fprintf('Done. ');
            toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % extract spikes and calculate firing rates
            firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
            firings_all{iter} = firings;
            bins = -warmup/1000:binSize:totalT/1000;
            [~,frLAct] = spikes2FR(firings,bins,inds.actLE); firingRates{iter,5} = frLAct;
            [~,frRAct] = spikes2FR(firings,bins,inds.actRE); firingRates{iter,6} = frRAct;

        end

        % threshold firing rates and score performance
        corrChoices = reshape(double(cellfun(@(x)ismember(x,{'Maltose','Octaacetate'}),stimuli)),[],1);
        firingRates_thresh = cell(numTrials,6);
        frLAct_max = max(cellfun(@max,firingRates(:,5)));
        frRAct_max = max(cellfun(@max,firingRates(:,6)));
        for TR = 1:numTrials    
            firingRates_thresh{TR,5} = threshold(firingRates{TR,5},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
            firingRates_thresh{TR,6} = threshold(firingRates{TR,6},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
        end

        choices = NaN(numTrials,1);
        bin1 = floor((stim.tStart+bound1)/binSize)+1; 
        bin2 = ceil((stim.tStart+bound2)/binSize);
        for trial = 1:numTrials
            if any(firingRates_thresh{trial,5}(bin1:bin2)) && ~any(firingRates_thresh{trial,6}(bin1:bin2))
                choices(trial) = 0;
            elseif ~any(firingRates_thresh{trial,5}(bin1:bin2)) && any(firingRates_thresh{trial,6}(bin1:bin2))
                choices(trial) = 1;
            end
        end

        perf{session_count,center_count} = [sum(choices==corrChoices),sum(~isnan(choices)),numTrials];

    end

end

%save(sprintf('perf_silStr%i_silDur%i_stimFall%i_stimGain%i.mat',...
%    round(silence_strength),round(100*silence_duration),round(stim_tauFall*1000),round(stim_Gain*100)),'perf');
    
%% PART 2

% options and global parameters -------------------------
sessions = 1:10;
numTrials = 100;

silence_center = 1.5; % [s]
silence_duration = 0.25; % [s]
silence_strength = 100; % [% increase]

warmup = 0; % simulate for this long before considering the trial to have begun, [ms]
totalT = 3000; % trial time, [ms]
dt = 0.05;

% main simulations -------------------------
timev = 0:dt:warmup+totalT;
windowPre = 1000*[silence_center-3*silence_duration/2,silence_center-silence_duration/2]+warmup;
windowDur = 1000*[silence_center-silence_duration/2,silence_center+silence_duration/2]+warmup;
allSpikes = cell(10,1);
frs = struct();

for session_count = sessions

    % session-specific parameters
    data = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,244+session_count));
    BinW = data.BinW;
    binSize = data.binSize;
    gate_startTimes = data.gate.StartTimes;
    gate_duration = data.gate.Duration;
    gate_floor = data.gate.Floor;
    gate_ceiling = data.gate.Ceiling;
    gate_switch = data.gate.Switch;
    neu = data.neu;
    net = data.net;
    stim = data.stim;
    inds = data.inds;
    S = data.S;
    clust = [data.clustSuc,data.clustQui,data.clustMal,data.clustOct,data.clustCueL,...
        data.clustCueR,data.clustActL,data.clustActR];
    clust = [clust,setdiff(1:14,clust)];
    
    silence_curve = getSilencing(silence_center,silence_duration,silence_strength,warmup,totalT,dt);

    fprintf('\n\nSession %i, silencing centered at %.3f:\n\n',244+session_count,silence_center);

    stimuli = [repmat({'Sucrose'},1,ceil(numTrials/4)),... 
               repmat({'Maltose'},1,ceil(numTrials/4)),...
               repmat({'Quinine'},1,ceil(numTrials/4)),...
               repmat({'Octaacetate'},1,ceil(numTrials/4))];
           
    preFRsI = zeros(14,numTrials); durFRsI = zeros(14,numTrials);
    preFRsE = zeros(14,numTrials); durFRsE = zeros(14,numTrials);

    for iter = 1:numTrials

        timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nSimulating trial %i...',iter);

        gate_curve = getActionGate(gate_startTimes(iter),gate_duration,gate_floor,gate_ceiling,...
            warmup,totalT,dt);
        [spkdata,~] = simulation(net,neu,stim,inds,stimuli{iter},gate_switch,gate_curve,S,true,...
            warmup,BinW,timev,silence_curve,false);

        fprintf('Done. ');
        toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extract spikes
        firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
        
        % calculate firing rates in windows
        for c = 1:length(clust)
            indI = net.indI(clust(c)):net.indI(clust(c)+1)-1;
            indE = net.indE(clust(c)):net.indE(clust(c)+1)-1;
            preFRsI(c,iter) = sum(ismember(firings(:,2),indI)&...
                firings(:,1)>=windowPre(1)&firings(:,1)<windowPre(2));
            durFRsI(c,iter) = sum(ismember(firings(:,2),indI)&...
                firings(:,1)>=windowDur(1)&firings(:,1)<windowDur(2));
            preFRsE(c,iter) = sum(ismember(firings(:,2),indE)&...
                firings(:,1)>=windowPre(1)&firings(:,1)<windowPre(2));
            durFRsE(c,iter) = sum(ismember(firings(:,2),indE)&...
                firings(:,1)>=windowDur(1)&firings(:,1)<windowDur(2));
        end      

    end

    frs(session_count).preI = preFRsI/71/silence_duration;
    frs(session_count).durI = durFRsI/71/silence_duration;
    frs(session_count).preE = preFRsE/250/silence_duration;
    frs(session_count).durE = durFRsE/250/silence_duration;

end

%save(sprintf('deltaFRs_str%i_dur%i_c%i',silence_strength,round(1000*silence_duration),round(1000*silence_center)),'frs');

%% PART 3

% options and global parameters -------------------------
session = 1;
trial = 1;

silence_center = 1.5;
silence_duration = 0.25; % [s]
silence_strength = 100; % [% increase]

binSize = 5; % [ms]
upsample_factor = 100;
smoothing = 800;

warmup = 0; % simulate for this long before considering the trial to have begun, [ms]
totalT = 3000; % trial time, [ms]
dt = 0.05;

% simulation -------------------------
rng(99995);
timev = 0:dt:warmup+totalT;
% session-specific parameters
data = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,244+session));
BinW = data.BinW;
inds = data.inds;
clustersWithRoles = data.clustersWithRoles;
gate_startTimes = data.gate.StartTimes;
gate_duration = data.gate.Duration;
gate_floor = data.gate.Floor;
gate_ceiling = data.gate.Ceiling;
gate_switch = data.gate.Switch;
net = data.net;
neu = data.neu;
stim = data.stim;
S = data.S;

silence_curve = getSilencing(silence_center,silence_duration,silence_strength,warmup,totalT,dt);

stimuli = [repmat({'Sucrose'},1,25),... 
           repmat({'Maltose'},1,25),...
           repmat({'Quinine'},1,25),...
           repmat({'Octaacetate'},1,25)];

gate_curve = getActionGate(gate_startTimes(trial),gate_duration,gate_floor,gate_ceiling,warmup,...
    totalT,dt);
[spkdata,~] = simulation(net,neu,stim,inds,stimuli{trial},gate_switch,gate_curve,S,true,warmup,...
    BinW,timev,silence_curve,false);

% extract spikes and calculate firing rates
firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
bins = -warmup/1000:binSize/1000:totalT/1000;
clust = [data.clustSuc,data.clustQui,data.clustMal,data.clustOct,data.clustCueL,data.clustCueR,...
    data.clustActL,data.clustActR];
clust = [clust,setdiff(1:14,clust)];
fr = struct();
for c = 1:length(clust)
    indE = net.indE(clust(c)):(net.indE(clust(c)+1)-1);
    indI = net.indI(clust(c)):(net.indI(clust(c)+1)-1);
    [t,fr(c).E] = spikes2FR(firings,bins,indE);
    [t,fr(c).I] = spikes2FR(firings,bins,indI);
end

t = myUpsample(t,upsample_factor);

figure(1); clf; hold all;
stim_curve = stim.Curve(1:length(timev));
plot(timev/1000,stim_curve/max(stim_curve)*100,'color',[0.5,0.5,0.5],'linewidth',1);
plot(timev/1000,gate_curve*100,'--','color',[0.5,0.5,0.5],'linewidth',1);
plot(timev/1000,100*(silence_curve-1)/max(silence_curve-1),...
    'color',[0.9290 0.6940 0.1250],'linewidth',1.5);
patch([silence_center-silence_duration/2,silence_center+silence_duration/2,...
    silence_center+silence_duration/2,silence_center-silence_duration/2],...
    [0,0,100,100],[0.9290 0.6940 0.1250],...
    'edgecolor','none','facecolor',[0.9290 0.6940 0.1250],'facealpha',0.5);
set(gca,'TickDir','out','color','none','box','off','fontsize',16);
xlim([0,3]);

figure(2); clf; hold all;
rasterPlotAbbrev(firings,inds,clustersWithRoles,totalT,net,7,true,'tick');
patch([silence_center-silence_duration/2,silence_center+silence_duration/2,...
    silence_center+silence_duration/2,silence_center-silence_duration/2],...
    [-5,-5,5005,5005],[0.9290 0.6940 0.1250],...
    'edgecolor','none','facecolor',[0.9290 0.6940 0.1250],'facealpha',0.4);
xlim([0,3]);
              
figure(3); clf; hold all;
for c = 1:length(clust)
    x = myUpsample(fr(c).E,upsample_factor);
    [xfilt,tfilt] = fun.gaussfilt(t,x,smoothing);
    if ismember(c,[5,7])
        color = 'g';
    elseif ismember(c,[6,8])
        color = 'b';
    else
        color = 'k';
    end
    plot(tfilt,xfilt,color,'linewidth',1.5);
end
myYLim = ylim();
patch([silence_center-silence_duration/2,silence_center+silence_duration/2,...
    silence_center+silence_duration/2,silence_center-silence_duration/2],...
    [myYLim(1),myYLim(1),myYLim(2),myYLim(2)],[0.9290 0.6940 0.1250],...
    'edgecolor','none','facecolor',[0.9290 0.6940 0.1250],'facealpha',0.3);
xlabel('Time [s]','interpreter','latex','fontsize',16);
ylabel('E cluster firing rates [Hz]','interpreter','latex','fontsize',16);
xlim([silence_center-1.5*silence_duration,silence_center+1.5*silence_duration]);
xticks(silence_center-1.5*silence_duration:silence_duration/2:silence_center+1.5*silence_duration);
ylim(myYLim);
set(gca,'TickDir','out','Color','none','box','off','fontsize',16);

figure(4); clf; hold all;
for c = 1:length(clust)
    x = myUpsample(fr(c).I,upsample_factor);
    [xfilt,tfilt] = fun.gaussfilt(t,x,smoothing);
    plot(tfilt,xfilt,'r','linewidth',1.5);
end
myYLim = ylim();
patch([silence_center-silence_duration/2,silence_center+silence_duration/2,...
    silence_center+silence_duration/2,silence_center-silence_duration/2],...
    [myYLim(1),myYLim(1),myYLim(2),myYLim(2)],[0.9290 0.6940 0.1250],...
    'edgecolor','none','facecolor',[0.9290 0.6940 0.1250],'facealpha',0.3);
xlabel('Time [s]','interpreter','latex','fontsize',16);
ylabel('I cluster firing rates [Hz]','interpreter','latex','fontsize',16);
xlim([silence_center-1.5*silence_duration,silence_center+1.5*silence_duration]);
xticks(silence_center-1.5*silence_duration:silence_duration/2:silence_center+1.5*silence_duration);
ylim(myYLim);
set(gca,'TickDir','out','Color','none','box','off','fontsize',16);

%% PART 4

% options and global parameters -------------------------
factorInOmissions = true;
isSmooth = true;
slopeType = 'endpoints0point5'; % 'midpoints','endpoints','endpoints0point5','regression'
upsample_factor = 100;
smoothing = 250;

% main -------------------------
% load baseline data
baseline_perf = cell(10,1);
i = 0;
for session = 245:254
    i = i + 1;
    data = load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,session));
    baseline_perf{i} = [sum(data.corrChoices==data.choices),sum(~isnan(data.choices)),100];
end
if factorInOmissions
    baseline_acc = cellfun(@(x)x(1)+0.5*(100-x(2)),baseline_perf);
else
    baseline_acc = cellfun(@(x)x(1)/x(2)*100,baseline_perf);
end
load(sprintf('%s/HMMData/Simulation/onsetTimeDist_QCS.mat',pwd));
load(sprintf('%s/HMMData/Simulation/onsetTimeDist_CCS.mat',pwd));
load(sprintf('%s/HMMData/Simulation/onsetTimeDist_ACS.mat',pwd));

% plot 1 (effect of varying silence strengths on performance): load
silence_strengths = [25,50,75,100,125]; % [%]
silence_center = 0:0.05:3;
curves = zeros(length(silence_strengths),length(silence_center));
for i = 1:length(silence_strengths)
    data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr%i_silDur250_stimFall160_stimGain60.mat',pwd,silence_strengths(i)));
    if factorInOmissions
        opto_acc = cellfun(@(x)x(1)+0.5*(100-x(2)),data.perf);
    else
        opto_acc = cellfun(@(x)x(1)/x(2)*100,data.perf);
    end
    curves(i,:) = mean(opto_acc);
end

% plot 1 (effect of varying silence strengths on performance): plot
figure(5); clf; hold all;
myYLim = [45,90];
colors = fun.distinguishable_colors(length(silence_strengths));
patch([3*onsetTimeDist_QCS(2,:),fliplr(3*onsetTimeDist_QCS(2,:))],...
    [myYLim(1)*ones(1,size(onsetTimeDist_QCS,2)),...
    fliplr(onsetTimeDist_QCS(1,:)*diff(myYLim)/max(onsetTimeDist_QCS(1,:))+myYLim(1))],...
    'r','EdgeColor','none','FaceAlpha',0.3);
patch([3*onsetTimeDist_CCS(2,:),fliplr(3*onsetTimeDist_CCS(2,:))],...
    [myYLim(1)*ones(1,size(onsetTimeDist_CCS,2)),...
    fliplr(onsetTimeDist_CCS(1,:)*diff(myYLim)/max(onsetTimeDist_CCS(1,:))+myYLim(1))],...
    'c','EdgeColor','none','FaceAlpha',0.3);
patch([3*onsetTimeDist_ACS(2,:),fliplr(3*onsetTimeDist_ACS(2,:))],...
    [myYLim(1)*ones(1,size(onsetTimeDist_ACS,2)),...
    fliplr(onsetTimeDist_ACS(1,:)*diff(myYLim)/max(onsetTimeDist_ACS(1,:))+myYLim(1))],...
    'b','EdgeColor','none','FaceAlpha',0.3);
patch([0,3,3,0],[mean(baseline_acc)-sqrt(var(baseline_acc)),mean(baseline_acc)-sqrt(var(baseline_acc)),...
    mean(baseline_acc)+sqrt(var(baseline_acc)),mean(baseline_acc)+sqrt(var(baseline_acc))],...
    [0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3);
h = [];
for i = 1:length(silence_strengths)
    x = silence_center; y = curves(i,:);
    if isSmooth
        x = myUpsample(x,upsample_factor); y = myUpsample(y,upsample_factor);
        [y,x] = fun.gaussfilt(x,y,smoothing);
    end
    h(i) = plot(x,y,'color',colors(i,:),'linewidth',2.5);
end
xlim([0,3]); ylim(myYLim);
set(gca,'TickDir','out','color','none','box','off','ticklabelinterpreter','latex','fontsize',14);
xlabel('Center of silencing window [s]','interpreter','latex','fontsize',18);
ylabel('Mean network accuracy [\%]','interpreter','latex','fontsize',18);
legend(h,arrayfun(@(x)['Strength = ',num2str(x),'\%'],silence_strengths,'UniformOutput',false),...
    'location','northeastoutside','interpreter','latex');

% plots 2,3,4 (effect of varying stimuli on performance): load
tau_Fs = [160,280,400,455,520,705];
silence_center = 0:0.05:3;
silence_center = silence_center(1:ceil(length(silence_center)/2));
curves = zeros(length(tau_Fs),length(silence_center));
for i = 1:length(tau_Fs)
    data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDur250_stimFall%i_stimGain60.mat',pwd,tau_Fs(i)));
    if factorInOmissions
        opto_acc = cellfun(@(x)x(1)+0.5*(100-x(2)),data.perf);
    else
        opto_acc = cellfun(@(x)x(1)/x(2)*100,data.perf);
    end
    curves(i,:) = mean(opto_acc(:,1:length(silence_center)));
end

% plot 2 (effect of varying stimuli on performance): plot
figure(6); clf; hold all;
colors = fun.distinguishable_colors(length(tau_Fs));
myYLim = [45,85];
h = [];
for i = 1:length(tau_Fs)
    x = silence_center; y = curves(i,:);
    if isSmooth
        x = myUpsample(x,upsample_factor); y = myUpsample(y,upsample_factor);
        [y,x] = fun.gaussfilt(x,y,smoothing);
    end
    h(i) = plot(x,y,'color',colors(i,:),'linewidth',2.5);
end
xlim([0,silence_center(end)]); ylim(myYLim);
set(gca,'TickDir','out','color','none','box','off','ticklabelinterpreter','latex','fontsize',14);
ylabel('Mean network accuracy [\%]','interpreter','latex','fontsize',18);
labels = arrayfun(@(x)['$\tau_\mathrm{F}=',num2str(x),'\ \mathrm{ms}$'],tau_Fs,'UniformOutput',false);
legend(h,labels,'interpreter','latex','fontsize',14,'location','northeastoutside');

% plot 3 (varying stimuli profiles): plot
figure(7); clf; hold all;
h = [];
t = 0:0.01:silence_center(end);
for i = 1:length(tau_Fs)
    s = (exp(-t/(tau_Fs(i)/1000))-exp(-t/0.150))/((tau_Fs(i)/1000)-0.150);
    s = s/max(s)*60;
    h(i) = plot(t,s,'linewidth',2.5,'color',colors(i,:));
end
xlim([0,silence_center(end)]); ylim([0,60]);
set(gca,'TickDir','out','color','none','box','off','ticklabelinterpreter','latex','fontsize',14);
ylabel('Stimulus gain [\%]','interpreter','latex','fontsize',18);
xlabel('Center of silencing window [s]','interpreter','latex','fontsize',18);
labels = arrayfun(@(x)['$\tau_\mathrm{F}=',num2str(x),'\ \mathrm{ms}$'],tau_Fs,'UniformOutput',false);
legend(h,labels,'interpreter','latex','fontsize',14,'location','northeastoutside');

% plot 4 (slope vs. tau): plot
figure(8); clf; hold all;
myXLim = [100,765];
slopes = [];
for i = 1:length(tau_Fs)
    if strcmp(slopeType,'midpoints')
        slope_i = NaN(1,length(curves(i,:))-2);
        for j = 2:length(curves(i,:))-1
            slope_i(j-1) = mean([(curves(i,j)-curves(i,j-1))/0.05,...
                (curves(i,j+1)-curves(i,j))/0.05]);
        end
        slope_i_mean = mean(slope_i);
    elseif strcmp(slopeType,'endpoints')
        slope_i_mean = (curves(i,end)-curves(i,1))/(0.05*(length(curves(i,:))-1));
    elseif strcmp(slopeType,'endpoints0point5')
        slope_i_mean = (curves(i,silence_center==0.5)-curves(i,1))/0.5;
    elseif strcmp(slopeType,'regression')
        params = polyfit(silence_center,curves(i,:),1);
        slope_i_mean = params(1);
    end
    slopes = [slopes,slope_i_mean];
    plot(tau_Fs(i),slope_i_mean,'.','markersize',20,'color',colors(i,:));
end
params = polyfit(tau_Fs(1:end),slopes,1);
plot(myXLim(1):0.01:myXLim(2),polyval(params,myXLim(1):0.01:myXLim(2)),':k','linewidth',2.5);
R2 = 1 - sum((slopes-polyval(params,tau_Fs(1:end))).^2)/sum((slopes-mean(slopes)).^2);
r2 = corr(tau_Fs(1:end)',slopes')^2;
set(gca,'TickDir','out','color','none','box','off','ticklabelinterpreter','latex','fontsize',14);
ylabel('Average accuracy slope [\%/s]','interpreter','latex','fontsize',18);
xlabel('$\tau_D$ [ms]','interpreter','latex','fontsize',18);
xlim(myXLim);
myYLim = ylim();
text(myXLim(1)+0.75*diff(myXLim),myYLim(1)+0.25*diff(myYLim),sprintf('$r^2=%.3f$',r2),'interpreter','latex');

% plot 5 (effect of varying silencing strengths on firing rates): load and plot
figure(9); clf; hold all;
silence_strengths = [25,50,75,100,125];
clusterAves = NaN(2,5);
for s = 1:length(silence_strengths)
    data = load(sprintf('%s/ProcessedData/Simulation/deltaFRs_str%i_dur250_c1500.mat',pwd,silence_strengths(s)));
    frs = data.frs;
    durFRsI = []; preFRsI = [];
    durFRsE = []; preFRsE = [];
    for n = 1:10
        durFRsI = [durFRsI,frs(n).durI];
        preFRsI = [preFRsI,frs(n).preI];
        durFRsE = [durFRsE,frs(n).durE];
        preFRsE = [preFRsE,frs(n).preE];
    end
    %yI = (mean(durFRsI,2));
    %yE = (mean(durFRsE,2));
    yI = (mean(durFRsI,2)-mean(preFRsI,2))./mean(preFRsI,2)*100;
    yE = (mean(durFRsE,2)-mean(preFRsE,2))./mean(preFRsE,2)*100;
    %yI(7:8) = NaN;
    plot(silence_strengths(s)+6*rand(1,14)-3,yI,'.r','markersize',15);
    plot(silence_strengths(s)+6*rand(1,14)-3,yE,'.k','markersize',15);
    plot(silence_strengths(s),median(yI(~isnan(yI)&~isinf(yI))),'or','markersize',15);
    plot(silence_strengths(s),median(yE(~isnan(yE)&~isinf(yE))),'ok','markersize',15);
    clusterAves(1,s) = median(yE(~isnan(yE)&~isinf(yE)));
    clusterAves(2,s) = median(yI(~isnan(yI)&~isinf(yI)));
end
yline(0,':k');
plot(silence_strengths,clusterAves(1,:),'--k','linewidth',1);
plot(silence_strengths,clusterAves(2,:),'--r','linewidth',1);
xticks(silence_strengths); xlim([20,130]);
set(gca,'tickdir','out','color','none','box','off');
xlabel('Silencing strength [\%]','interpreter','latex','fontsize',14);
ylabel('$\Delta f_c^{rel}$ [\%]','interpreter','latex','fontsize',14);
title(['$\Delta f_c^{rel}=\frac{\mathrm{mean}_{n,i}\big(f_{c,n,i}^{post}\big)-',...
    '\mathrm{mean}_{n,i}\big(f_{c,n,i}^{pre}\big)}{\mathrm{mean}_{n,i}',...
    '\big(f_{c,n,i}^{pre}\big)}\times100\%$'],'interpreter','latex','fontsize',18);

%% local function definitions

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

function [SPIKE_DATA,INPUT_DATA] = simulation(NET,NEU,STIM,INDS,STIMULUS,GATE_SWITCH,...
    GATE_CURVE,S,IS_MODIFY_WEIGHTS,WARMUP,BIN_W,TIME_V,SILENCING,SAVE_INPUTS)
% initialization
q = 1; fired = []; binfired = [];
V = 0 + 4*randn(NET.N,1); % initial membrane potential
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

function rasterPlotAbbrev(SPIKE_DATA,INDS,CLUSTERS_WITH_ROLES,TOTAL_T,NET,DOWNSAMPLE,APPROX_RATIO,DOT_OR_TICK)
downsampleEC = DOWNSAMPLE;
anchor = 0;
if APPROX_RATIO
    downsampleIC = max(round(downsampleEC/250*71),1);
    downsampleEB = 2*downsampleEC;
else
    downsampleIC = downsampleEC;
    downsampleEB = downsampleEC;
end 
% plot E taste clusters
ind_net = [INDS.sucE(randperm(length(INDS.sucE),downsampleEC)),...
    INDS.quiE(randperm(length(INDS.quiE),downsampleEC)),...
    INDS.malE(randperm(length(INDS.malE),downsampleEC)),...
    INDS.octE(randperm(length(INDS.octE),downsampleEC))];
ind_remap = (1:4*downsampleEC)'+anchor;
anchor = anchor + 4*downsampleEC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'k.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'k');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2,'Suc','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+downsampleEC,'Qui','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+2*downsampleEC,'Mal','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+3*downsampleEC,'Oct','VerticalAlignment','middle','fontsize',12);
% plot E cue left and action left clusters
ind_net = [INDS.cueLE(randperm(length(INDS.cueLE),downsampleEC)),...
    INDS.actLE(randperm(length(INDS.actLE),downsampleEC))];
ind_remap = [(1:downsampleEC)';((2*downsampleEC+1):3*downsampleEC)']+anchor;
anchor = anchor + downsampleEC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'g.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'g');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+4*downsampleEC,'Cue L','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+6*downsampleEC,'Act L','VerticalAlignment','middle','fontsize',12);
% plot E cue right and action right clusters
ind_net = [INDS.cueRE(randperm(length(INDS.cueRE),downsampleEC)),...
    INDS.actRE(randperm(length(INDS.actRE),downsampleEC))];
ind_remap = [(1:downsampleEC)';((2*downsampleEC+1):3*downsampleEC)']+anchor;
anchor = anchor + 3*downsampleEC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'b.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'b');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+5*downsampleEC,'Cue R','VerticalAlignment','middle','fontsize',12);
text(TOTAL_T/1000+0.05,(1+downsampleEC)/2+7*downsampleEC,'Act R','VerticalAlignment','middle','fontsize',12);
% plot remaining E clusters
for clust = setdiff(1:NET.Q,CLUSTERS_WITH_ROLES)
    ind_net = NET.indE(clust):(NET.indE(clust+1)-1); 
    ind_net = ind_net(randperm(length(ind_net),downsampleEC));
    ind_remap = (1:downsampleEC)'+anchor;
    anchor = anchor + downsampleEC;
    times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
    ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
    if strcmp(DOT_OR_TICK,'dot')
        plot(times,ind_plot,'k.','markersize',0.1);
    elseif strcmp(DOT_OR_TICK,'tick')
        times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
        plot([times;times],[ind_plot-0.4;ind_plot+0.4],'k');
    else
        error('Invalid option for DOT_OR_TICK. Use lowercase letters');
    end
end
% plot background E neurons
ind_net = NET.indE(NET.Q+1):(NET.indE(NET.Q+2)-1);
ind_net = ind_net(randperm(length(ind_net),downsampleEB));
ind_remap = (1:downsampleEB)'+anchor;
anchor = anchor + downsampleEB;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'k.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'k');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot I taste clusters
ind_net = [INDS.sucI(randperm(length(INDS.sucI),downsampleIC)),...
    INDS.quiI(randperm(length(INDS.quiI),downsampleIC)),...
    INDS.malI(randperm(length(INDS.malI),downsampleIC)),...
    INDS.octI(randperm(length(INDS.octI),downsampleIC))];
ind_remap = (1:4*downsampleIC)'+anchor;
anchor = anchor + 4*downsampleIC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'r.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot I cue left and action left clusters
ind_net = [INDS.cueLI(randperm(length(INDS.cueLI),downsampleIC)),...
    INDS.actLI(randperm(length(INDS.actLI),downsampleIC))];
ind_remap = [(1:downsampleIC)';((2*downsampleIC+1):3*downsampleIC)']+anchor;
anchor = anchor + downsampleIC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'r.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot I cue right and action right clusters
ind_net = [INDS.cueRI(randperm(length(INDS.cueRI),downsampleIC)),...
    INDS.actRI(randperm(length(INDS.actRI),downsampleIC))];
ind_remap = [(1:downsampleIC)';((2*downsampleIC+1):3*downsampleIC)']+anchor;
anchor = anchor + 3*downsampleIC;
times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
if strcmp(DOT_OR_TICK,'dot')
    plot(times,ind_plot,'r.','markersize',0.1);
elseif strcmp(DOT_OR_TICK,'tick')
    times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
    plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
else
    error('Invalid option for DOT_OR_TICK. Use lowercase letters');
end
% plot remaining I clusters
for clust = setdiff(1:NET.Q,CLUSTERS_WITH_ROLES)
    ind_net = NET.indI(clust):(NET.indI(clust+1)-1); 
    ind_net = ind_net(randperm(length(ind_net),downsampleIC));
    ind_remap = (1:downsampleIC)'+anchor;
    anchor = anchor + downsampleIC;
    times = SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),1)/1000;
    ind_plot = arrayfun(@(x)ind_remap(x==ind_net),SPIKE_DATA(ismember(SPIKE_DATA(:,2),ind_net),2));
    if strcmp(DOT_OR_TICK,'dot')
        plot(times,ind_plot,'r.','markersize',0.1);
    elseif strcmp(DOT_OR_TICK,'tick')
        times = reshape(times,1,[]); ind_plot = reshape(ind_plot,1,[]);
        plot([times;times],[ind_plot-0.4;ind_plot+0.4],'r');
    else
        error('Invalid option for DOT_OR_TICK. Use lowercase letters');
    end
end
% visual formatting
yline(NET.Q*downsampleEC+downsampleEB+0.5,'r');
yline(NET.Q*(downsampleEC+downsampleIC)+downsampleEB+0.5,'r');
set(gca,'fontsize',18,'tickdir','out','color','none','box','off');
yticks([]);
xlim([0, TOTAL_T/1000]); ylim([0,NET.Q*(downsampleEC+downsampleIC)+downsampleEB+1]);
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

function SILENCING_CURVE = getSilencing(CENTER,DURATION,STRENGTH,WARMUP,TOTAL_T,DT)
timeV = 0:DT:(WARMUP+TOTAL_T);
start = CENTER-DURATION/2;
SILENCING_CURVE = ones(1,length(timeV));
SILENCING_CURVE((timeV>=WARMUP+1000*start)&(timeV<=WARMUP+1000*start+1000*DURATION)) = 1 + STRENGTH/100;
end

function UPSAMPLED_VEC = myUpsample(VEC,NUM)
UPSAMPLED_VEC = [];
for i = 1:(length(VEC)-1)
    UPSAMPLED_VEC = [UPSAMPLED_VEC,linspace(VEC(i),VEC(i+1),NUM+2)];
end
end