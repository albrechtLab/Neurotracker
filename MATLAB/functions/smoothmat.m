%% smoothmat
%   apply smoothing across rows or columns

function output = smoothmat(M,dim,span)
% USAGE: function output = smoothmat(M,dim,span)
%
%       operates the smooth function on rows or columns of matrix M
%       according to dimension dim.  span default is 5.

    if nargin < 3 span = 5; end
    if nargin < 2 dim = 1; end

    sm = size(M);

    if dim > 1 M = M'; sm = size(M); end

    smoothM = [];
    for i = 1:sm(1)
        smoothM = [smoothM; smooth(M(i,:),span)'];
    end

    if dim > 1
        output = smoothM'; 
    else
        output = smoothM;
    end
end