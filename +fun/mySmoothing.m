function [X_DATA_SMOOTH,Y_DATA_SMOOTH] = mySmoothing(X_DATA,Y_DATA,PADDING_TYPE,NUM_PADS,EXPANSION,FILTER_CONSTANT)
%[X_DATA_SMOOTH,Y_DATA_SMOOTH] = mySmoothing(X_DATA,Y_DATA,PADDING_TYPE,NUM_PADS,EXPANSION,FILTER_CONSTANT)
%
%Smooth a function for plotting
%
%    Input X_DATA: vector, the 'x-axis' data of your function
%
%    Input Y_DATA: vector, the 'y-axis' data of your function
%
%    Input PADDING_TYPE: string 
%        'zeros': 0s will be appended to the left and right sides of your x
%                 and y data before smoothing
%        'continuation': x data will be continued linearly on both sides, y
%                        data will have its first value appended to the
%                        left and its last value appended to the right
%
%    Input NUM_PADS: integer, the number of pad values to append on each
%    side of the data (if this argument is 0, there will be no padding,
%    regardless of the previous argument)
%
%    Input EXPANSION: integer, point density after linear interpolation
%    (for any two adjacent points in the original data, EXPANSION is the
%    number of points between them in the smoothed data, including
%    themselves --- therefore an EXPANSION value of 2 or less means no
%    linear interpolation is performed)
%
%    Input FILTER_CONSTANT: integer, binwidth of the Gaussian kernel that
%    gaussfilt uses for smoothing (e.g. 10 means the Gaussian's width
%    stretches over 10 consecutive x values)
%
% -LL
%

% padding
x = [fliplr(X_DATA(1)-(X_DATA(2)-X_DATA(1))*(1:NUM_PADS)),X_DATA,X_DATA(end)+(X_DATA(end)-X_DATA(end-1))*(1:NUM_PADS)];
if strcmp(PADDING_TYPE,'zeros')
    y = [zeros(1,NUM_PADS),Y_DATA,zeros(1,NUM_PADS)];
elseif strcmp(PADDING_TYPE,'continuation')
    y = [Y_DATA(1)*ones(1,NUM_PADS),Y_DATA,Y_DATA(end)*ones(1,NUM_PADS)];
end
if EXPANSION > 1
    % linear interpolation
    x_interp = x(1); y_interp = y(1);
    for i = 1:length(x)-1
        interpX = linspace(x(i),x(i+1),EXPANSION);
        interpY = linspace(y(i),y(i+1),EXPANSION);
        x_interp = [x_interp,interpX(2:end)];
        y_interp = [y_interp,interpY(2:end)];
    end
else
    x_interp = x;
    y_interp = y;
end
% gaussian filtering
[y_filt,x_filt] = fun.gaussfilt(x_interp,y_interp,FILTER_CONSTANT,0);
% return
X_DATA_SMOOTH = x_filt;
Y_DATA_SMOOTH = y_filt;
end