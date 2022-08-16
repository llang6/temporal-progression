%% plotStateTimes
%
% Generate the distribution of onset times plot for experiment or
% simulation data (Figures 2 and 5).
%
% - LL
%

%% dependencies
% requires access to files:
% /HMMData/Experiment/classifiedStates_exp.mat
% /HMMData/Experiment/classifiedStates_exp_shuff_circ.mat
% /HMMData/Experiment/classifiedStates_exp_shuff_swap.mat
% /HMMData/Simulation/classifiedStates_sim.mat
% /HMMData/Simulation/classifiedStates_sim_shuff_circ.mat
% /HMMData/Simulation/classifiedStates_sim_shuff_swap.mat
% requires access to functions:
% loadVar (in +fun)
% hist2curve (in +fun)
% gaussfilt (in +fun)

%% parameters
expOrSim = 'simulation'; % 'experiment', 'simulation'
dataType = 'original'; % 'original', 'shuffled_circular', 'shuffled_swap'
excludedSessions = [5,7,17,19,20]; % experiment sessions 5,7,17,19,20 were excluded from analysis
adjustForPaddingInClassification = true;
clipOutsideRange = true;
correctTrialsOnly = false;
filterConstant = 400; 
showHistograms = false;
binWidth = 0.05; % histogram bin width (units of warped time)
numPads = 2; % number of 0 bins to add to each side of histogram
expansion = 100; % linear interpolation factor for connecting histogram peaks

%% setup, load data, and format for plotting
homeDir = pwd; addpath(homeDir);
model = sprintf('%s_%s',expOrSim,dataType);
if strcmp(expOrSim,'simulation'), excludedSessions = []; end
cd('HMMData');
switch model
    case 'simulation_original'
        cd('Simulation'); classifiedStates = fun.loadVar('classifiedStates_sim.mat');
    case 'simulation_shuffled_circular'
        cd('Simulation'); classifiedStates = fun.loadVar('classifiedStates_sim_shuff_circ.mat');
    case 'simulation_shuffled_swap'
        cd('Simulation'); classifiedStates = fun.loadVar('classifiedStates_sim_shuff_swap.mat');
    case 'experiment_original'
        cd('Experiment'); classifiedStates = fun.loadVar('classifiedStates_exp.mat');
    case 'experiment_shuffled_circular'
        cd('Experiment'); classifiedStates = fun.loadVar('classifiedStates_exp_shuff_circ.mat');
    case 'experiment_shuffled_swap'
        cd('Experiment'); classifiedStates = fun.loadVar('classifiedStates_exp_shuff_swap.mat');
end
cd(homeDir);
QCS_onsets = []; QCS_offsets = []; QCS_count = 0;
CCS_onsets = []; CCS_offsets = []; CCS_count = 0;
ACS_onsets = []; ACS_offsets = []; ACS_count = 0;
for i = 2:size(classifiedStates,1)
    if ismember(classifiedStates{i,1},excludedSessions), continue; end
    if strcmp(classifiedStates{i,3},'Exclusive Quality-coding')
        if correctTrialsOnly
            ind = ismember(classifiedStates{i,8},find(classifiedStates{i,9}==1));
        else
            ind = true(1,length(classifiedStates{i,8}));
        end
        QCS_onsets = [QCS_onsets,classifiedStates{i,6}(ind)];
        QCS_offsets = [QCS_offsets,classifiedStates{i,7}(ind)];
        QCS_count = QCS_count + 1;
    elseif strcmp(classifiedStates{i,3},'Cue-coding')
        if correctTrialsOnly
            ind = ismember(classifiedStates{i,8},find(classifiedStates{i,9}==1));
        else
            ind = true(1,length(classifiedStates{i,8}));
        end
        CCS_onsets = [CCS_onsets,classifiedStates{i,6}(ind)];
        CCS_offsets = [CCS_offsets,classifiedStates{i,7}(ind)];
        CCS_count = CCS_count + 1;
    elseif strcmp(classifiedStates{i,3},'Action-coding')
        if correctTrialsOnly
            ind = ismember(classifiedStates{i,8},find(classifiedStates{i,9}==1));
        else
            ind = true(1,length(classifiedStates{i,8}));
        end
        ACS_onsets = [ACS_onsets,classifiedStates{i,6}(ind)];
        ACS_offsets = [ACS_offsets,classifiedStates{i,7}(ind)];
        ACS_count = ACS_count + 1;
    end
end
if contains(model,'experiment')
    if adjustForPaddingInClassification
        flags = find(((QCS_onsets<-1)&(QCS_offsets<-1))|((QCS_onsets>0)&(QCS_offsets>0)));
        QCS_onsets(flags) = []; QCS_offsets(flags) = [];
        flags = find(((CCS_onsets<-1)&(CCS_offsets<-1))|((CCS_onsets>0)&(CCS_offsets>0)));
        CCS_onsets(flags) = []; CCS_offsets(flags) = [];
        flags = find(((ACS_onsets<-1)&(ACS_offsets<-1))|((ACS_onsets>0)&(ACS_offsets>0)));
        ACS_onsets(flags) = []; ACS_offsets(flags) = [];
    end
    if clipOutsideRange
        QCS_onsets_adj = max(min(QCS_onsets,0),-1)+1;
        CCS_onsets_adj = max(min(CCS_onsets,0),-1)+1;
        ACS_onsets_adj = max(min(ACS_onsets,0),-1)+1;
        QCS_offsets_adj = max(min(QCS_offsets,0),-1)+1;
        CCS_offsets_adj = max(min(CCS_offsets,0),-1)+1;
        ACS_offsets_adj = max(min(ACS_offsets,0),-1)+1;
    else
        QCS_onsets_adj = QCS_onsets+1;
        CCS_onsets_adj = CCS_onsets+1;
        ACS_onsets_adj = ACS_onsets+1;
        QCS_offsets_adj = QCS_offsets+1;
        CCS_offsets_adj = CCS_offsets+1;
        ACS_offsets_adj = ACS_offsets+1;
    end
else
    QCS_onsets_adj = QCS_onsets;
    CCS_onsets_adj = CCS_onsets;
    ACS_onsets_adj = ACS_onsets;
    QCS_offsets_adj = QCS_offsets;
    CCS_offsets_adj = CCS_offsets;
    ACS_offsets_adj = ACS_offsets;
end

%% plot distributions of onset times
figureDimensions = [689,239]; 

f = figure(1); clf;
f.Position(3:4) = figureDimensions;
bins = -0.1:binWidth:1.1;
% first distribution
h1 = histogram(QCS_onsets_adj,bins,'normalization','pdf','EdgeColor','none','FaceColor','r'); hold on;
c1 = fun.hist2curve(h1,'numPads',numPads,'linewidth',2,'color','r','linestyle',':');
xData = c1.XData; yData = c1.YData;
if expansion > 0
    xData_interp = xData(1); yData_interp = yData(1);
    for i = 1:length(xData)-1
        interpX = linspace(xData(i),xData(i+1),expansion);
        interpY = linspace(yData(i),yData(i+1),expansion);
        xData_interp = [xData_interp,interpX(2:end)];
        yData_interp = [yData_interp,interpY(2:end)];
    end    
else
    xData_interp = xData; yData_interp = yData;
end
[c1filt,t1filt] = fun.gaussfilt(xData_interp,yData_interp,filterConstant,0);
% onsetTimeDist_QCS = [c1filt' ; t1filt]; save('onsetTimeDist_QCS.mat','onsetTimeDist_QCS');
p1 = plot(t1filt,c1filt,'linewidth',2.5,'color','r'); hold on;
if ~showHistograms, h1.Visible = 'off'; end
c1.Visible = 'off';
% second distribution
h2 = histogram(CCS_onsets_adj,bins,'normalization','pdf','EdgeColor','none','FaceColor','c'); hold on;
c2 = fun.hist2curve(h2,'numPads',numPads,'linewidth',2,'color','c','linestyle',':'); hold on;
xData = c2.XData; yData = c2.YData;
if expansion > 0
    xData_interp = xData(1); yData_interp = yData(1);
    for i = 1:length(xData)-1
        interpX = linspace(xData(i),xData(i+1),expansion);
        interpY = linspace(yData(i),yData(i+1),expansion);
        xData_interp = [xData_interp,interpX(2:end)];
        yData_interp = [yData_interp,interpY(2:end)];
    end    
else
    xData_interp = xData; yData_interp = yData;
end
[c2filt,t2filt] = fun.gaussfilt(xData_interp,yData_interp,filterConstant,0);
% onsetTimeDist_CCS = [c2filt' ; t2filt]; save('onsetTimeDist_CCS.mat','onsetTimeDist_CCS');
p2 = plot(t2filt,c2filt,'linewidth',2.5,'color','c'); hold on; 
if ~showHistograms, h2.Visible = 'off'; end 
c2.Visible = 'off';
% third distribution
h3 = histogram(ACS_onsets_adj,bins,'normalization','pdf','EdgeColor','none','FaceColor','b'); hold on;
c3 = fun.hist2curve(h3,'numPads',numPads,'linewidth',2,'color','b','linestyle',':'); hold on;
xData = c3.XData; yData = c3.YData;
if expansion > 0
    xData_interp = xData(1); yData_interp = yData(1);
    for i = 1:length(xData)-1
        interpX = linspace(xData(i),xData(i+1),expansion);
        interpY = linspace(yData(i),yData(i+1),expansion);
        xData_interp = [xData_interp,interpX(2:end)];
        yData_interp = [yData_interp,interpY(2:end)];
    end    
else
    xData_interp = xData; yData_interp = yData;
end
[c3filt,t3filt] = fun.gaussfilt(xData_interp,yData_interp,filterConstant,0);
% onsetTimeDist_ACS = [c3filt' ; t3filt]; save('onsetTimeDist_ACS.mat','onsetTimeDist_ACS');
p3 = plot(t3filt,c3filt,'linewidth',2.5,'color','b'); hold on; 
if ~showHistograms, h3.Visible = 'off'; end
c3.Visible = 'off';
legend([p1,p2,p3],{sprintf('Quality-coding states (N=%i from %i states)',length(QCS_onsets_adj),QCS_count),...
                   sprintf('Cue-coding states (N=%i from %i states)',length(CCS_onsets_adj),CCS_count),...
                   sprintf('Action-coding states (N=%i from %i states)',length(ACS_onsets_adj),ACS_count)},...
                   'location','eastoutside');
if contains(model,'simulation')
    xlabel('Time (relative units: 0 = Trial start, 1 = Trial end)'); 
else
    xlabel('Time (relative units: 0 = Taste, 1 = Decision)');
end
ylabel('Probability density');
xlim([-0.1,1.1]);
ylim([0,max([c1filt',c2filt',c3filt'])+0.1]);
temp = split(model,'_');
titleString = temp{1}; for i = 2:length(temp), titleString = [titleString,'\_',temp{i}]; end
title(sprintf('Distribution of onset times: %s',titleString));
set(gca,'fontsize',10,'Color','none','TickDir','out','box','off');
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'Renderer','painters');
% print(sprintf('onsetTimes_%s',model),'-depsc');
% print(sprintf('onsetTimes_%s',model),'-dpng');

fprintf('\nMeans\n-------\n');
fprintf('Quality-coding: %.2f\n',mean(QCS_onsets_adj));
fprintf('Cue-coding: %.2f\n',mean(CCS_onsets_adj));
fprintf('Action-coding: %.2f\n',mean(ACS_onsets_adj));
fprintf('\nInferred peaks\n-------\n');
fprintf('Quality-coding: %.2f\n',t1filt(c1filt==max(c1filt)));
fprintf('Cue-coding: %.2f\n',t2filt(c2filt==max(c2filt)));
fprintf('Action-coding: %.2f\n',t3filt(c3filt==max(c3filt)));