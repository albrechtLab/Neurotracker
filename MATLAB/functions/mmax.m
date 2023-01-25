%% mmax
%   multidimensional maximum

function output = mmax(input,dimlist)

    if nargin < 2 dimlist = 1:ndims(input); end

    temp = input;
    temp(isnan(temp)) = -Inf;

    for i = 1:length(dimlist)
        dimswap = 1:ndims(input); dimswap(dimlist(i)) = 1; dimswap(1) = dimlist(i);
        temp = permute(temp,dimswap);
        %temp = sum(temp,dimlist(i));
        temp = max(temp);
        temp = permute(temp,dimswap);
    end
    output = squeeze(temp);
end