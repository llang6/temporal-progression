%% plotSilencing_allNetworks
%
% Plot the data created with optoSilencing_allNetworks
%
% This was used to create panels for Figure 6 and Supplementary Figures
% 5-6.
%
% -LL
%

%% dependencies
% requires access to data:
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain60_sameV0.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain60.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall705_stimGain60_sameV0.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall705_stimGain60.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0_removedTail.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_removedTail.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0_removedHead.mat
% /ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_removedHead.mat
% /ProcessedData/Simulation/perf_silStr100_silDur250_stimFall160_stimGain200.mat
% /ProcessedData/Simulation/perf_silStr25_silDur250_stimFall160_stimGain200.mat
% /RawData/Simulation/simulationData252_noSilencing_0_sameV0.mat

%% stimulus: 60%, 160 ms ("original")
% === load ===
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain60_sameV0.mat',pwd));
temp = vertcat(data.perf{:,1});
perf_baseline_V0 = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain60.mat',pwd));
temp = vertcat(data.perf{:,2});
perf_silencing_samp = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
temp = vertcat(data.perf{:,3});
perf_silencing_delay = temp(:,1)+0.5*(temp(:,3)-temp(:,2));

% === stats (1-factor within-subjects ANOVA) ===
% setup according to: http://compneurosci.com/wiki/images/e/e6/Repeated_ANOVA_MATLAB_v2.pdf
% perf is (number of subjects)-by-(number of conditions)
% column order:
% full stim (gain: 60%, tau: 160 ms), no silencing;
% full stim (gain: 60%, tau: 160 ms), silencing during samp;
% full stim (gain: 60%, tau: 160 ms), silencing during delay;
perf = NaN(10,3);
perf(:,1) = perf_baseline_V0;
perf(:,2) = perf_silencing_samp;
perf(:,3) = perf_silencing_delay;
% name the variables by column number
varNames = cell(size(perf,2),1);
for i = 1:size(perf,2)
    varNames{i} = strcat('V',num2str(i));
end
% create a table storing the responses
tPerf = array2table(perf,'VariableNames',varNames);
% create a table reflecting the within-subject factors
silencing = cell(size(perf,2),1); % silencing conditions
% assign values to the parameters based on the data sorting
silencing{1} = 'none';
silencing{2} = 'sampling';
silencing{3} = 'delay';
% create the within table
factorNames = {'Silencing'};
within = table(silencing,'VariableNames',factorNames);
% fit the repeated measures model
rm = fitrm(tPerf,sprintf('V1-V%i~1',size(perf,2)),'WithinDesign',within);
ranovatblb = ranova(rm,'WithinModel','Silencing');
% post-hoc test
posthoc_bonferroni = multcompare(rm,'Silencing','ComparisonType','bonferroni');

% === plots ===
% --- parameters ---
fontsize_small = 10;
fontsize_big = 12;
width = 2.25;
height = 1.94;

figure(1); clf;
% --- plot performance --- 
subplot(100,3,121:300); hold all;
plot([ones(1,length(perf_baseline_V0));2*ones(1,length(perf_baseline_V0));3*ones(1,length(perf_baseline_V0))],...
    [perf_baseline_V0';perf_silencing_samp';perf_silencing_delay'],'k','linewidth',0.5);
plot(1,perf_baseline_V0,'.k','markersize',10);
plot(2,perf_silencing_samp,'.r','markersize',10);
plot(3,perf_silencing_delay,'.r','markersize',10);
plot(1,mean(perf_baseline_V0),'*k','markersize',10);
plot(2,mean(perf_silencing_samp),'*r','markersize',10);
plot(3,mean(perf_silencing_delay),'*r','markersize',10);
patch([0,4,4,0],[mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),...
    mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0))],...
    [0.4,0.4,0.4],'edgecolor','none','facealpha',0.3);
xlim([0.7,3.3]);
xticks(1:3);
xticklabels({'No silencing','Silencing: [0, 0.5]','Silencing: [0.5, 3]'});
ylim([40,100]);
yticks(40:10:100);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Network performance [%]','fontsize',fontsize_big);

% --- plot inputs ---
timev = 0:0.05:3000;
tau_D = 160;
Gain = 60;
stim = (1/(tau_D-150))*(exp(-timev/tau_D) - exp(-timev/150));
stim = stim/max(abs(stim))*Gain; 

subplot(100,3,1:3:58); hold all;
plot(timev/1000,stim,'k','linewidth',1);
yticks([0,30,60]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Stimulus gain [%]','fontsize',fontsize_big);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,2:3:59); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0,0.5,0.5,0],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,30,60]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,3:3:60); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0.5,3,3,0.5],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
yticks([0,30,60]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

set(gcf,'Renderer','painters'); 
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[width,height]); 
set(gcf,'PaperPosition',[0,0,width,height]); 
set(gcf,'Units','inches'); 

%print(gcf,'stimulus_original_silencing_asInExp','-depsc');

%% stimulus: 60%, 705 ms ("longer")
% === load ===
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall705_stimGain60_sameV0.mat',pwd));
temp = vertcat(data.perf{:,1});
perf_baseline_V0 = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall705_stimGain60.mat',pwd));
temp = vertcat(data.perf{:,2});
perf_silencing_samp = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
temp = vertcat(data.perf{:,3});
perf_silencing_delay = temp(:,1)+0.5*(temp(:,3)-temp(:,2));

% === stats ===
perf = NaN(10,3);
perf(:,1) = perf_baseline_V0;
perf(:,2) = perf_silencing_samp;
perf(:,3) = perf_silencing_delay;
varNames = cell(size(perf,2),1);
for i = 1:size(perf,2)
    varNames{i} = strcat('V',num2str(i));
end
tPerf = array2table(perf,'VariableNames',varNames);
silencing = cell(size(perf,2),1); 
silencing{1} = 'none';
silencing{2} = 'sampling';
silencing{3} = 'delay';
factorNames = {'Silencing'};
within = table(silencing,'VariableNames',factorNames);
rm = fitrm(tPerf,sprintf('V1-V%i~1',size(perf,2)),'WithinDesign',within);
ranovatblb = ranova(rm,'WithinModel','Silencing');
posthoc_bonferroni = multcompare(rm,'Silencing','ComparisonType','bonferroni');

% === plots ===
% --- parameters ---
fontsize_small = 10;
fontsize_big = 12;
width = 2.25;
height = 1.94;

figure(2); clf;
% --- plot performance --- 
subplot(100,3,121:300); hold all;
plot([ones(1,length(perf_baseline_V0));2*ones(1,length(perf_baseline_V0));3*ones(1,length(perf_baseline_V0))],...
    [perf_baseline_V0';perf_silencing_samp';perf_silencing_delay'],'k','linewidth',0.5);
plot(1,perf_baseline_V0,'.k','markersize',10);
plot(2,perf_silencing_samp,'.r','markersize',10);
plot(3,perf_silencing_delay,'.r','markersize',10);
plot(1,mean(perf_baseline_V0),'*k','markersize',10);
plot(2,mean(perf_silencing_samp),'*r','markersize',10);
plot(3,mean(perf_silencing_delay),'*r','markersize',10);
patch([0,4,4,0],[mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),...
    mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0))],...
    [0.4,0.4,0.4],'edgecolor','none','facealpha',0.3);
xlim([0.7,3.3]);
xticks(1:3);
xticklabels({'No silencing','Silencing: [0, 0.5]','Silencing: [0.5, 3]'});
ylim([40,100]);
yticks(40:10:100);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Network performance [%]','fontsize',fontsize_big);

% --- plot inputs ---
timev = 0:0.05:3000;
tau_D = 705;
Gain = 60;
stim = (1/(tau_D-150))*(exp(-timev/tau_D) - exp(-timev/150));
stim = stim/max(abs(stim))*Gain; 

subplot(100,3,1:3:58); hold all;
plot(timev/1000,stim,'k','linewidth',1);
yticks([0,30,60]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Stimulus gain [%]','fontsize',fontsize_big);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,2:3:59); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0,0.5,0.5,0],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,30,60]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,3:3:60); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0.5,3,3,0.5],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,30,60]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

set(gcf,'Renderer','painters'); 
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[width,height]); 
set(gcf,'PaperPosition',[0,0,width,height]); 
set(gcf,'Units','inches');  
 
%print(gcf,'stimulus_longer705_silencing_asInExp','-depsc');

%% stimulus: 200%, 160 ms ("stronger")
% === load ===
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0.mat',pwd));
temp = vertcat(data.perf{:,1});
perf_baseline_V0 = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200.mat',pwd));
temp = vertcat(data.perf{:,2});
perf_silencing_samp = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
temp = vertcat(data.perf{:,3});
perf_silencing_delay = temp(:,1)+0.5*(temp(:,3)-temp(:,2));

% === plots ===
% --- parameters ---
fontsize_small = 10;
fontsize_big = 12;
width = 2.25;
height = 1.94;

figure(3); clf;
% --- plot performance --- 
subplot(100,3,121:300); hold all;
plot([ones(1,length(perf_baseline_V0));2*ones(1,length(perf_baseline_V0));3*ones(1,length(perf_baseline_V0))],...
    [perf_baseline_V0';perf_silencing_samp';perf_silencing_delay'],'k','linewidth',0.5);
plot(1,perf_baseline_V0,'.k','markersize',10);
plot(2,perf_silencing_samp,'.r','markersize',10);
plot(3,perf_silencing_delay,'.r','markersize',10);
plot(1,mean(perf_baseline_V0),'*k','markersize',10);
plot(2,mean(perf_silencing_samp),'*r','markersize',10);
plot(3,mean(perf_silencing_delay),'*r','markersize',10);
patch([0,4,4,0],[mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),...
    mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0))],...
    [0.4,0.4,0.4],'edgecolor','none','facealpha',0.3);
xlim([0.7,3.3]);
xticks(1:3);
xticklabels({'No silencing','Silencing: [0, 0.5]','Silencing: [0, 0.5]'});
ylim([40,100]);
yticks(40:10:100);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Network performance [%]','fontsize',fontsize_big);

% --- plot inputs ---
timev = 0:0.05:3000;
tau_D = 160;
Gain = 200;
stim = (1/(tau_D-150))*(exp(-timev/tau_D) - exp(-timev/150));
stim = stim/max(abs(stim))*Gain; 

subplot(100,3,1:3:58); hold all;
plot(timev/1000,stim,'k','linewidth',1);
yticks([0,100,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Stimulus gain [%]','fontsize',fontsize_big);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,2:3:59); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0,0.5,0.5,0],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,100,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,3:3:60); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0.5,3,3,0.5],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,100,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

set(gcf,'Renderer','painters'); 
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[width, height]);  
set(gcf,'PaperPosition',[0,0,width,height]); 
set(gcf,'Units','inches'); 

% print(gcf,'stimulus_stronger200_silencing_asInExp','-depsc');

%% stimulus: 200%, 160 ms, head only
% === load ===
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0_removedTail.mat',pwd));
temp = vertcat(data.perf{:,1});
perf_baseline_V0 = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_removedTail.mat',pwd));
temp = vertcat(data.perf{:,2});
perf_silencing_samp = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
temp = vertcat(data.perf{:,3});
perf_silencing_delay = temp(:,1)+0.5*(temp(:,3)-temp(:,2));

% === plots ===
% --- parameters ---
fontsize_small = 10;
fontsize_big = 12;
width = 2.25;
height = 1.94;

figure(4); clf;
% --- plot performance --- 
subplot(100,3,121:300); hold all;
plot([ones(1,length(perf_baseline_V0));2*ones(1,length(perf_baseline_V0));3*ones(1,length(perf_baseline_V0))],...
    [perf_baseline_V0';perf_silencing_samp';perf_silencing_delay'],'k','linewidth',0.5);
plot(1,perf_baseline_V0,'.k','markersize',10);
plot(2,perf_silencing_samp,'.r','markersize',10);
plot(3,perf_silencing_delay,'.r','markersize',10);
plot(1,mean(perf_baseline_V0),'*k','markersize',10);
plot(2,mean(perf_silencing_samp),'*r','markersize',10);
plot(3,mean(perf_silencing_delay),'*r','markersize',10);
patch([0,4,4,0],[mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),...
    mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0))],...
    [0.4,0.4,0.4],'edgecolor','none','facealpha',0.3);
xlim([0.7,3.3]);
xticks(1:3);
xticklabels({'No silencing','Silencing: [0, 0.5]','Silencing: [0.5, 3]'});
ylim([40,100]);
yticks(40:10:100);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Network performance [%]','fontsize',fontsize_big);

% --- plot inputs ---
timev = 0:0.05:3000;
tau_D = 160;
Gain = 200;
stim = (1/(tau_D-150))*(exp(-timev/tau_D) - exp(-timev/150));
stim = stim/max(abs(stim))*Gain; 
stim(timev>=500) = 0;

subplot(100,3,1:3:58); hold all;
plot(timev/1000,stim,'k','linewidth',1);
yticks([0,100,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Stimulus gain [%]','fontsize',fontsize_big);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,2:3:59); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0,0.5,0.5,0],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,100,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,3:3:60); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0.5,3,3,0.5],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,100,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

set(gcf,'Renderer','painters'); 
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[width,height]);  
set(gcf,'PaperPosition',[0,0,width,height]); 
set(gcf,'Units','inches'); 

%print(gcf,'stimulus_stronger200_removedTail_silencing_asInExp','-depsc');

%% stimulus: 200%, 160 ms, tail only
% === load ===
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0_removedHead.mat',pwd));
temp = vertcat(data.perf{:,1});
perf_baseline_V0 = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_removedHead.mat',pwd));
temp = vertcat(data.perf{:,2});
perf_silencing_samp = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
temp = vertcat(data.perf{:,3});
perf_silencing_delay = temp(:,1)+0.5*(temp(:,3)-temp(:,2));

% === plots ===
% --- parameters ---
fontsize_small = 10;
fontsize_big = 12;
width = 2.25;
height = 1.94;

figure(5); clf;
% --- plot performance --- 
subplot(100,3,121:300); hold all;
plot([ones(1,length(perf_baseline_V0));2*ones(1,length(perf_baseline_V0));3*ones(1,length(perf_baseline_V0))],...
    [perf_baseline_V0';perf_silencing_samp';perf_silencing_delay'],'k','linewidth',0.5);
plot(1,perf_baseline_V0,'.k','markersize',10);
plot(2,perf_silencing_samp,'.r','markersize',10);
plot(3,perf_silencing_delay,'.r','markersize',10);
plot(1,mean(perf_baseline_V0),'*k','markersize',10);
plot(2,mean(perf_silencing_samp),'*r','markersize',10);
plot(3,mean(perf_silencing_delay),'*r','markersize',10);
patch([0,4,4,0],[mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)-sqrt(var(perf_baseline_V0)),...
    mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0)),mean(perf_baseline_V0)+sqrt(var(perf_baseline_V0))],...
    [0.4,0.4,0.4],'edgecolor','none','facealpha',0.3);
xlim([0.7,3.3]);
xticks(1:3);
xticklabels({'No silencing';'Silencing: [0,0.5]';'Silencing: [0.5,3]'});
ylim([40,100]);
yticks(40:10:100);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Network performance [%]','fontsize',fontsize_big);

% --- plot inputs ---
timev = 0:0.05:3000;
tau_D = 160;
Gain = 200;
stim = (1/(tau_D-150))*(exp(-timev/tau_D) - exp(-timev/150));
stim = stim/max(abs(stim))*Gain; 
stim(timev<=500) = 0;

subplot(100,3,1:3:58); hold all;
plot(timev/1000,stim,'k','linewidth',1);
yticks([0,100,200]);
ylim([0,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
ylabel('Stimulus gain [%]','fontsize',fontsize_big);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,2:3:59); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0,0.5,0.5,0],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,100,200]);
ylim([0,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

subplot(100,3,3:3:60); hold all;
plot(timev/1000,stim,'k','linewidth',1);
patch([0.5,3,3,0.5],[0,0,Gain,Gain],'y','edgecolor','none','facealpha',0.4);
%yline(stim(timev==500),':k');
yticks([0,100,200]);
ylim([0,200]);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Time [s]','fontsize',fontsize_big);

set(gcf,'Renderer','painters'); 
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[width,height]); 
set(gcf,'PaperPosition',[0,0,width,height]); 
set(gcf,'Units','inches'); 

%print(gcf,'stimulus_stronger200_removedHead_silencing_asInExp','-depsc');

%% 2-factor within-subjects ANOVA for stimulus: 200%, 160 ms 
% perf is (number of subjects)-by-(number of conditions)
% column order:
% full stim, no silencing;
% full stim, silencing during samp;
% full stim, silencing during delay;
% head stim, no silencing;
% head stim, silencing during samp;
% head stim, silencing during delay;
% tail stim, no silencing;
% tail stim, silencing during samp;
% tail stim, silencing during delay

perf = NaN(10,9);
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0.mat',pwd));
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,1));
perf(:,1) = data_column;
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200.mat',pwd));
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,2));
perf(:,2) = data_column;
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,3));
perf(:,3) = data_column;
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0_removedTail.mat',pwd));
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,1));
perf(:,4) = data_column;
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_removedTail.mat',pwd));
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,2));
perf(:,5) = data_column;
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,3));
perf(:,6) = data_column;
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0_removedHead.mat',pwd));
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,1));
perf(:,7) = data_column;
data = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_removedHead.mat',pwd));
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,2));
perf(:,8) = data_column;
data_column = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data.perf(:,3));
perf(:,9) = data_column;

varNames = cell(size(perf,2),1);
for i = 1:size(perf,2)
    varNames{i} = strcat('V',num2str(i));
end

% create a table storing the responses
tPerf = array2table(perf,'VariableNames',varNames);

% create a table reflecting the within-subject factors
stimulus = cell(size(perf,2),1); % stimulus conditions
silencing = cell(size(perf,2),1); % silencing conditions

% assign values to the parameters based on the data sorting
stimulus(1:3) = repmat({'full'},3,1);
stimulus(4:6) = repmat({'head_only'},3,1);
stimulus(7:9) = repmat({'tail_only'},3,1);
silencing([1,4,7]) = repmat({'none'},3,1);
silencing([2,5,8]) = repmat({'sampling'},3,1);
silencing([3,6,9]) = repmat({'delay'},3,1);

% create the within table
factorNames = {'Stimulus','Silencing'};
within = table(stimulus,silencing,'VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(tPerf,sprintf('V1-V%i~1',size(perf,2)),'WithinDesign',within);
ranovatblb = ranova(rm,'WithinModel','Stimulus*Silencing');

posthoc1_bonferroni = multcompare(rm,'Stimulus','By','Silencing','ComparisonType','bonferroni');
posthoc2_bonferroni = multcompare(rm,'Silencing','By','Stimulus','ComparisonType','bonferroni');

%% moving window silencing (of strength 25% or 100%) for networks with stimulus: 200%, 160 ms
% === load ===
data_baseline = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDurExp_stimFall160_stimGain200_sameV0.mat',pwd)); % known V0
temp = vertcat(data_baseline.perf{:,1});
perf_baseline = temp(:,1)+0.5*(temp(:,3)-temp(:,2));
% [~,ind] = min(abs(perf_baseline-mean(perf_baseline))); % ind = 8;
data_100 = load(sprintf('%s/ProcessedData/Simulation/perf_silStr100_silDur250_stimFall160_stimGain200.mat',pwd));
acc_100 = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data_100.perf);
data_25 = load(sprintf('%s/ProcessedData/Simulation/perf_silStr25_silDur250_stimFall160_stimGain200.mat',pwd));
acc_25 = cellfun(@(x)x(1)+0.5*(x(3)-x(2)),data_25.perf);

% === plot ===
% --- parameters ---
linewidth = 1;
markersize = 10;
fontsize_vsmall = 8;
fontsize_small = 10;
fontsize_big = 12;
width = 6.8;
height = 4.2;

figure(6); clf; hold all;
p1 = plot((0:60)/60*3,mean(acc_100,1),'r','linewidth',linewidth);
plot((0:60)/60*3,mean(acc_100,1),'.r','markersize',markersize);
p2 = plot((0:60)/60*3,mean(acc_25,1),'k','linewidth',linewidth);
plot((0:60)/60*3,mean(acc_25,1),'.k','markersize',markersize);
% mean +/- 1 stdev for control condition
patch([0,3,3,0],[mean(perf_baseline)-sqrt(var(perf_baseline)),mean(perf_baseline)-sqrt(var(perf_baseline)),...
    mean(perf_baseline)+sqrt(var(perf_baseline)),mean(perf_baseline)+sqrt(var(perf_baseline))],[0.4,0.4,0.4],...
    'edgecolor','none','facealpha',0.4);
% "beginning"
patch([0,0.25,0.25,0],[45,45,90,90],'y','facealpha',0.4);
% "cue onset" (for single network 252)
data_raw_noSilencing = load(sprintf('%s/RawData/Simulation/simulationData252_noSilencing_0_sameV0.mat',pwd));
cueOnsetTimes = data_raw_noSilencing.cueOnsetTimes;
patch([mean(cueOnsetTimes)-0.125,mean(cueOnsetTimes)+0.125,mean(cueOnsetTimes)+0.125,mean(cueOnsetTimes)-0.125],...
    [45,45,90,90],'y','facealpha',0.4);
% "middle"
patch([1.375,1.625,1.625,1.375],[45,45,90,90],'y','facealpha',0.4);
xlim([0,3]);
ylim([45,90]);
yticks(45:5:90);
legend([p1,p2],{'100%','25%'},'location','northeast','fontsize',fontsize_vsmall);
set(gca,'tickdir','out','color','none','box','off','fontsize',fontsize_small);
xlabel('Center of silencing window [s]','fontsize',fontsize_big);
ylabel('Mean network accuracy [%]','fontsize',fontsize_big);

set(gcf,'Renderer','painters'); 
set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[width,height]); 
set(gcf,'PaperPosition',[0,0,width,height]); 
set(gcf,'Units','inches'); 

% print(gcf,'accuracyVsCenter_strongStim','-depsc');
