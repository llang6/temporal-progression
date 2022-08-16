function varargout = hist2curve(HISTOGRAM,varargin)
%hist2curve(HISTOGRAM)
%
% Takes a HISTOGRAM object and plots the curve connecting the peaks of the
% bin heights.
% PH = hist2curve(HISTOGRAM) returns the plot handle for the curve. 
%
% Optional Name,Value pair arguments:
% 'NumPads' - number of 0s to place on both sides of HISTOGRAM x values
% 'LineWidth' - passed to 'plot()'
% 'Color' - passed to 'plot()'
% 'LineStyle' - passed to 'plot()' 
%
% Remember to use 'hold on' between each histogram you create, or the
% previous one may be destroyed.
% Also remember not to close the plot window or the histogram object may be
% destroyed.
%
% Updated 11/30/21 to include LineStyle and NumPads as optional arguments
% and require all optional arguments to be given as name, value pairs.
%
% -LL
%

% parse input
if isempty(varargin)
    % default parameters
    NumPads = 2; 
    LineWidth = 2; 
    Color = 'k'; 
    LineStyle = '-';
else
    ind = find(cellfun(@(x)strcmpi(x,'NumPads'),varargin),1);
    if ~isempty(ind)
        NumPads = varargin{ind+1};
    else
        NumPads = 2;
    end
    ind = find(cellfun(@(x)strcmpi(x,'LineWidth'),varargin),1);
    if ~isempty(ind)
        LineWidth = varargin{ind+1};
    else
        LineWidth = 2;
    end
    ind = find(cellfun(@(x)strcmpi(x,'Color'),varargin),1);
    if ~isempty(ind)
        Color = varargin{ind+1};
    else
        Color = 'k';
    end
    ind = find(cellfun(@(x)strcmpi(x,'LineStyle'),varargin),1);
    if ~isempty(ind)
        LineStyle = varargin{ind+1};
    else
        LineStyle = '-';
    end
end

% main
x = HISTOGRAM.BinEdges(1:HISTOGRAM.NumBins) + HISTOGRAM.BinWidth/2;
y = HISTOGRAM.Values;
% Pad with zeros outside the support
x1 = []; x2 = [];
if NumPads > 0
    for i = 1:NumPads
        x1 = [x(1)-i*HISTOGRAM.BinWidth,x1];
        x2 = [x2,x(end)+i*HISTOGRAM.BinWidth];
    end
end
x = [x1,x,x2];
y = [zeros(1,NumPads),y,zeros(1,NumPads)];
% x = [x(1)-2*histogram.BinWidth, x(1)-histogram.BinWidth, ...
%      x, ...
%      x(end)+histogram.BinWidth, x(end)+2*histogram.BinWidth];
% y = [0, 0, y, 0, 0];
% To do: smooth?
ph = plot(x,y,'LineWidth',LineWidth,'Color',Color,'LineStyle',LineStyle);
if nargout > 0, varargout{1} = ph; end

end