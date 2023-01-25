%% msum
%   multidimensional nansum

function output = msum(input,dimlist)

    if nargin < 2 dimlist = 1:ndims(input); end %#ok<*SEPEX>

    temp = input;

    for i = 1:length(dimlist)
        dimswap = 1:ndims(input); dimswap(dimlist(i)) = 1; dimswap(1) = dimlist(i);
        temp = permute(temp,dimswap);
        %temp = sum(temp,dimlist(i));
        if size(temp,1) > 1
            temp = nansum(temp);
        end
        temp = permute(temp,dimswap);
    end
    output = squeeze(temp);
end
