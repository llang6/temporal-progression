%% optoSilencing_allNetworks
%
% This script is similar to lifNet_altStim, insofar as it loads the
% original 10 networks and runs new simulations using a modified form of
% the stimulus. However, it also introduces optogenetic silencing during
% the simulations. It is potentially very simulation-intensive, and only
% performance data is saved (spiking data is not). The loop over sessions
% was run as a parfor loop for speed, and this requires formatting
% variables in some particular ways. Here the parfor loop is commented out
% and replaced with a for loop for demonstration.
%
% Setting 'sil_type' to 'asInExp' will simulate under 3 conditions: no
% silencing, silencing during sampling [0,0.5], and silencing during delay
% [0.5,3]. This was used along with changing 'stim_tauFall' and 'stim_Gain'
% to find a stimulus resulting in performance results that matched the
% experiment (Figure 6 and Supplementary Figure 5).
%
% Setting 'sil_type' to 'movingWindow' will simulate under 61 conditions,
% with a 250 ms-wide silencing window centered in 50 ms increments from 0
% to 3 s. Used for Supplementary Figure 6. 
%
% You can choose to use the same initial conditions for all neurons in all
% trials across all conditions by setting 'conserveInitCond' to true. This
% will load the initial conditions from a static file 'V0.mat'. Ultimately
% for 'asInExp' simulations we used the same initial conditions for
% control/no silencing conditions across different stimuli, but not for
% different silencing conditions. For 'movingWindow' simulations we did not
% use the same initial conditions.
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Simulation/V0.mat
%
% requires access to functions:
% loadVar (in +fun)
% getSilencing (in +fun)
% getActionGate (in +fun)
% simulation (in +fun)
% rasterPlotAbbrev (in +fun) 
% spikes2FR (in +fun)

%% parameters
stim_tauFall = 0.160; % [s]
stim_Gain = 2.00; % [*100% increase]
stim_removeHead = 0;
stim_removeTail = 0;
sil_strength = 100;
sil_type = 'asInExp'; % 'asInExp', 'movingWindow' (also 'sampling' and 'delay' for plots)
conserveInitCond = true;
isPlot = true;

%% simulations
% setup
addpath(pwd); % gain access to package functions

% other fixed parameters 
sessions = 1:10;
numTrials = 100;
warmup = 0; % simulate for this long before considering the trial to have begun, [ms]
totalT = 3000; % trial time, [ms]
dt = 0.05; % [ms]
stim_tauRise = 0.150; % [s]

switch sil_type
    case 'asInExp'
        silence_center = [0,0.25,1.75]; % relative to stimulus onset, [s]
        silence_duration = [0,0.5,2.5]; % [s]
        silence_strength = [0,100,100]; % [% increase]
        append0 = 'Exp';
    case 'sampling'
        silence_center = 0.25; 
        silence_duration = 0.5; 
        silence_strength = 100; 
        append0 = 'Sampling';
    case 'delay'
        silence_center = 1.75; 
        silence_duration = 2.5;
        silence_strength = 100; 
        append0 = 'Delay';
    case 'movingWindow'
        silence_center = 0:0.05:3; % relative to stimulus onset, [s]
        silence_duration = 0.25*ones(1,length(silence_center)); % [s]
        silence_strength = sil_strength*ones(1,length(silence_center)); % [% increase]
        append0 = '250';
end

% main simulations -------------------------
timev = 0:dt:warmup+totalT;
index2 = 1:length(silence_center);
stim_Curve(timev>=warmup) = (1/((stim_tauFall*1000)-(stim_tauRise*1000)))*...
    (exp(-(0:dt:totalT)/(stim_tauFall*1000))-...
    exp(-(0:dt:totalT)/(stim_tauRise*1000)));
stim_Curve = stim_Curve/max(abs(stim_Curve))*stim_Gain;
if stim_removeTail
    stim_Curve(timev>=warmup+500) = 0;
end
if stim_removeHead
    stim_Curve(timev<warmup+500) = 0;
end
if conserveInitCond
    V0 = fun.loadVar(sprintf('%s/RawData/Simulation/V0.mat',pwd));
end

if isPlot
    figure(1); clf;
    plot(timev/1000,stim_Curve*100,'k','linewidth',2);
    set(gca,'tickdir','out','color','none','box','off');
end

perf = cell(10,length(silence_center));

%parfor session_count = sessions
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
    stim.Curve = stim_Curve;
    inds = data.inds;
    S = data.S;
    stateThresh = data.stateThresh;
    stateMinTime = data.stateMinTime;
    bound1 = data.bound1;
    bound2 = data.bound2;

    for center_count = index2
        
        silence_curve = fun.getSilencing(silence_center(center_count),silence_duration(center_count),...
            silence_strength(center_count),warmup,totalT,dt,'center');

        fprintf('\n\nSession %i, silencing (%i%%) from %.3f to %.3f...\n\n',244+session_count,...
            silence_strength(center_count),silence_center(center_count)-silence_duration(center_count)/2,...
            silence_center(center_count)+silence_duration(center_count)/2);

        firings_all = cell(numTrials,1);
        firingRates = cell(numTrials,6);
        stimuli = [repmat({'Sucrose'},1,ceil(numTrials/4)),... 
                   repmat({'Maltose'},1,ceil(numTrials/4)),...
                   repmat({'Quinine'},1,ceil(numTrials/4)),...
                   repmat({'Octaacetate'},1,ceil(numTrials/4))];

        for iter = 1:numTrials

            timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\nSimulating trial %i...',iter);

            gate_curve = fun.getActionGate(gate_startTimes(iter),gate_duration,gate_floor,gate_ceiling,warmup,totalT,dt);
            if conserveInitCond
                initCond = V0(session_count).V0(:,iter);
            else
                initCond = 4*randn(4994,1);
            end
            [spkdata,~] = fun.simulation(net,neu,stim,inds,stimuli{iter},gate_switch,gate_curve,S,true,...
                warmup,BinW,timev,silence_curve,false,initCond);

            fprintf('Done. ');
            toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % extract spikes and calculate firing rates
            firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
            firings_all{iter} = firings;
            if isPlot
                figure(2); clf; hold all;
                fun.rasterPlotAbbrev(firings,inds,data.clustersWithRoles,totalT,net,7,true,'tick',0.05);
                set(gca,'fontsize',12);
                if silence_duration(center_count)>0
                    patch([silence_center(center_count)-silence_duration(center_count)/2,...
                        silence_center(center_count)+silence_duration(center_count)/2,...
                        silence_center(center_count)+silence_duration(center_count)/2,...
                        silence_center(center_count)-silence_duration(center_count)/2],...
                        [-10,-10,net.N+1,net.N+1],'y','edgecolor','none','facealpha',0.4);
                end
                drawnow; pause;
            end
            bins = -warmup/1000:binSize:totalT/1000;
            [~,frLAct] = fun.spikes2FR(firings,bins,inds.actLE); firingRates{iter,5} = frLAct;
            [~,frRAct] = fun.spikes2FR(firings,bins,inds.actRE); firingRates{iter,6} = frRAct;

        end

        % threshold firing rates and score performance
        corrChoices = reshape(double(cellfun(@(x)ismember(x,{'Maltose','Octaacetate'}),stimuli)),[],1);
        firingRates_thresh = cell(numTrials,6);
        frLAct_max = max(cellfun(@max,firingRates(:,5)));
        frRAct_max = max(cellfun(@max,firingRates(:,6)));
        for TR = 1:numTrials  
            firingRates_thresh{TR,5} = fun.threshold(firingRates{TR,5},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
            firingRates_thresh{TR,6} = fun.threshold(firingRates{TR,6},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
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

% save performance data
if conserveInitCond, append1 = '_sameV0'; else, append1 = ''; end
if stim_removeHead, append2 = '_removedHead'; else, append2 = ''; end
if stim_removeTail, append3 = '_removedTail'; else, append3 = ''; end

% save(sprintf('%s/ProcessedData/Simulation/perf_silStr%i_silDur%s_stimFall%i_stimGain%i%s%s%s.mat',...
%    pwd,sil_strength,append0,round(stim_tauFall*1000),round(stim_Gain*100),append1,append2,append3),'perf');    
