% function [z tout] = gaussfilt(t,x,binwidth,causal_filter)
%
% filter 1D vector x with a Gaussian filter of std = binwidth (half width
% is half of that) by convolution of x with a gaussian. It uses 'conv'.
% Based on
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/272556
%
% INPUT:
% - t: time vector
% - x: data
% - binwidth: (optional) twice the halfwidth (=SD) of the gaussian kernel
% used for filtering, in units of time steps. NOTE: must be integer > 1;
% set it to 1 to prevent filtering. Default: 10 (time steps)
% - causal_filter: (flag; optional): if ==1 causal filter used (this means
% the filter=0 at negative times). Default: 0.
%
% OUTPUT: 
% - z: the filtered quantity
% - tout: appropriate time vector for z
%
% NOTES: a fictitious continuation of the t and x values to the left and
% right of the data vector are added to prevent bad filtering at the
% boundaries. 
% If no data are given, a dataset is generated for illustration. The data
% are meant to mimic the membrane potential of a spiking neuron. Change the
% value of 'causal_filter' inside code to visualize the difference. 
%
% GLC, Nov 7, 2013; revisions:
% - Feb 2017 (input arguments 'binwidth' and 'causal_filter' made optional
% and inverted the order of the outputs)
% - Jun 2020: added surrogate dataset if no data given.
% LL revisions:
% - Dec 2021: fixed output time vectors, which contained a repeated value
% that caused a 'kink' in the resulting plot, and which were improperly
% aligned for filters with even binwidths

function  [z, tout] = gaussfilt(t,x,varargin)

% parse input
isplot=0; % default
binwidth=10; % default: no filtering
if nargin>2 binwidth=varargin{1}; end
causal_filter=0; % default: a-causal filtering
if nargin>3 causal_filter=varargin{2}; end

% if no data given, generate some data for illustration
if nargin<2    
    t=0:0.1:10; % time in ms
    n=length(t); 
    x=randn(n,1); % example signal (gaussian 'membrane potential')
    x(randi(20,1,n)==1)=50; % 'spikes'; comment for a more standard example
    causal_filter=0; 
    isplot=1; % plot results for illustration
end

% preprocessing
if isrow(t) t=t'; end
if isrow(x) x=x'; end
exlen=binwidth;
x1=repmat(x(1),exlen,1);
x2=repmat(x(end),exlen,1);
x0=[x1 ; x ; x2];
tbin=diff(t(1:2));
t1=fliplr(t(1):-tbin:t(1)-(exlen-1)*tbin);
t2=fliplr(t(end):tbin:t(end)+(exlen-1)*tbin);
t0=[t1' ; t ; t2'];

% Construct blurring window:
windowWidth = binwidth;
halfWidth = ceil(windowWidth / 2);
gaussFilter = gausswin(binwidth);
if causal_filter gaussFilter(1:floor(max(size(gaussFilter))/2))=0; end
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

% Do the blur:
z = conv(x0, gaussFilter); 
z = z((halfWidth+exlen):(length(z)-halfWidth-exlen+1));
% tout = t0(halfWidth+exlen:end-halfWidth-exlen+1);
% tout=exlen*tbin+t0(1:length(z))-tbin;
% --- note on tout --------------------------------------------------------
% we want z(t) to be x(t) when the center of the filter was aligned with it
% but for filters with even binwidths, there is no center
% so, we take the data from half a binwidth before and after the full sweep
% this results in length(z) = length(x) + 1
% to realign t with z, add another bin and shift it back by half a bin
% -------------------------------------------------------------------------
if floor(binwidth/2)~=binwidth/2
    % if filter has odd binwidth length, use original time vector
    tout = t;
elseif floor(binwidth/2)==binwidth/2
    % if filter has even binwidth length, realign t to z
    tout = [t;t(end)+tbin]-tbin/2;
end

% plot
if isplot
    figure(101); clf;
    hold on;
    plot(t,x,'k','linewidth',2);
    plot(tout,z, 'r-', 'linewidth', 3);
    hold off;
    xlabel('time (ms)');
    ylabel('signal (a.u.)');
    set(gca,'fontsize',22,'color','none','fontweight','normal');
end
