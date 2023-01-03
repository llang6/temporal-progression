function rasterPlotAbbrev(SPIKE_DATA,INDS,CLUSTERS_WITH_ROLES,TOTAL_T,NET,DOWNSAMPLE,APPROX_RATIO,DOT_OR_TICK,X_OFFSET)
%rasterPlotAbbrev(SPIKE_DATA,INDS,CLUSTERS_WITH_ROLES,TOTAL_T,NET,DOWNSAMPLE,APPROX_RATIO,DOT_OR_TICK,X_OFFSET)
%
%Generate a raster plot with optional downsampling
%
%Make sure to grab a figure, clear it, and hold all before calling this
%function
% ex. figure(1); clf; hold all;
%     rasterPlotAbbrev(...);
%
%This function is very similar to @rasterPlot, though it offers greater
%flexibility. See that documentation for information about shared inputs
%SPIKE_DATA, INDS, CLUSTERS_WITH_ROLES, TOTAL_T, and NET.
%It does not require the INDEX_REMAP_INV argument that @rasterPlot does,
%and instead it does the index remapping manually here.
%
%Inputs unique to @rasterPlotAbbrev:
%    DOWNSAMPLE: integer, the number of neurons to plot for each E cluster
%
%    APPROX_RATIO: logical, if true the function will plot a number of
%        neurons for each I cluster and for the E background based on
%        DOWNSAMPLE and in proportion with their true numbers in the
%        population, if false the function will just plot DOWNSAMPLE
%        neurons for each entity
%
%    DOT_OR_TICK: string, if 'dot' spikes will be plotted as dots, if
%        'tick' spikes will be plotted as tick marks
%
%    X_OFFSET: number, controls how far from the right edge of the plot the
%        functional cluster labels wil be
%
% -LL
%

% calculate numbers of neurons to plot
downsampleEC = DOWNSAMPLE;
anchor = 0;
if APPROX_RATIO
    downsampleIC = max(round(downsampleEC/NET.Necl(1)*NET.Nicl(1)),1);
    downsampleEB = max(round(NET.Necl(end)/NET.Necl(1)*downsampleEC),1);
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
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2,'Suc','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+downsampleEC,'Qui','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+2*downsampleEC,'Mal','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+3*downsampleEC,'Oct','VerticalAlignment','middle','fontsize',10);
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
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+4*downsampleEC,'Cue L','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+6*downsampleEC,'Act L','VerticalAlignment','middle','fontsize',10);
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
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+5*downsampleEC,'Cue R','VerticalAlignment','middle','fontsize',10);
text(TOTAL_T/1000+X_OFFSET,(1+downsampleEC)/2+7*downsampleEC,'Act R','VerticalAlignment','middle','fontsize',10);
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
