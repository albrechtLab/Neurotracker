%% realTimeData 
%   reshape time-series data from stacked trials to real-time plots

function output = realTimeData(mat,exp,animal,interval)
% output = realTimeData(mat,exp,animal,interval)
%           mat = m x n matrix of time series data, m time points x n rows (trials) 
%           exp = n vector of experiment number
%           animal = n vector of animal number
%           interval = time points in real-time interval

    timepoints = size(mat,1);
    trials = size(mat,2);
    if nargin < 4 interval = timepoints; end
    
    if (interval < timepoints) error('Interval less than timepoints per trial'); end
    if (length(exp) ~= trials) error('experiment vector must be same length as matrix rows'); end
    if (length(animal) ~= trials) error('animal vector must be same length as matrix rows'); end

    alist = sort(unique(animal));
    elist = sort(unique(exp));
     
    RT.mat = NaN * zeros(interval*length(elist),length(alist));  % preallocate full matrix
    % loop through and build matrix row by row
    for r = 1:trials
        RT.mat((exp(r)-1)*interval + [1:timepoints],find(alist == animal(r))) = mat(:,r);
    end
    
    RT.animal = alist;              % new realtime animal rows
    RT.exp = ones(size(alist));     % all same experiment now
    
    output = RT;
end