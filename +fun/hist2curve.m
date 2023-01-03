function varargout = hist2curve(HISTOGRAM,varargin)
%hist2curve(HISTOGRAM)
%
% Takes a HISTOGRAM object and plots the curve connecting the peaks of the
% bin heights.
% PH = hist2curve(HISTOGRAM) returns the plot handle for the curve. 
% [PH,X,Y] = hist2curve(HISTOGRAM) returns the plot handle for the curve as
% well as the x and y values of the curve.
%
% Optional Name,Value pair arguments:
% 'NumPads' - number of 0s to place on both sides of HISTOGRAM x values
% 'isPlot' - whether or not you want the function to plot the curve
% 'LineWidth' - passed to 'plot()'
% 'Color' - passed to 'plot()'
% 'LineStyle' - passed to 'plot()' 
%
% Remember to use 'hold on' between each histogram you create, or the
% previous one may be destroyed.
% Also remember not to close the plot window or the histogram object may be
% destroyed.
% 
% Example use:
% myData = randn(1,100);
% figure(1); clf; hold all;
% myHist = histogram(myData);
% [ph,x,y] = hist2curve(myHist);
%
% Updated 12/29/22 to include isPlot as an optional argument and allow up
% to 3 returned values.
% 
% Updated 11/30/21 to include LineStyle and NumPads as optional arguments
% and require all optional arguments to be given as name, value pairs.
%
% -LL
%

% default parameters
default_NumPads = 2;
default_isPlot = true;
default_LineWidth = 2;
default_Color = 'k';
default_LineStyle = '-';

% parse input
if isempty(varargin)
    % use defaults
    NumPads = default_NumPads; 
    isPlot = default_isPlot;
    LineWidth = default_LineWidth; 
    Color = default_Color; 
    LineStyle = default_LineStyle;
else
    ind = find(cellfun(@(x)strcmpi(x,'NumPads'),varargin),1);
    if ~isempty(ind), NumPads = varargin{ind+1}; else, NumPads = default_NumPads; end
    ind = find(cellfun(@(x)strcmpi(x,'isPlot'),varargin),1);
    if ~isempty(ind), isPlot = varargin{ind+1}; else, isPlot = default_isPlot; end
    ind = find(cellfun(@(x)strcmpi(x,'LineWidth'),varargin),1);
    if ~isempty(ind), LineWidth = varargin{ind+1}; else, LineWidth = default_LineWidth; end
    ind = find(cellfun(@(x)strcmpi(x,'Color'),varargin),1);
    if ~isempty(ind), Color = varargin{ind+1}; else, Color = default_Color; end
    ind = find(cellfun(@(x)strcmpi(x,'LineStyle'),varargin),1);
    if ~isempty(ind), LineStyle = varargin{ind+1}; else, LineStyle = default_LineStyle; end
end

% --- main ---
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
% plot
if isPlot
    ph = plot(x,y,'LineWidth',LineWidth,'Color',Color,'LineStyle',LineStyle);
else
    ph = [];
end

% return
if nargout > 0
    varargout{1} = ph; 
    varargout{2} = x;
    varargout{3} = y;
end

end