function M2 = binaverage(M, x, crop)
%  function M2 = binaverage(M, x, crop)
%
%  where M is the matrix with row-data to be averaged
%    and x is the multiplier (# of bins to combine)
%    and crop = 1 drops any extra originalbins

    if nargin < 3 crop = 0; end

    M2 = [];
    lenM = size(M,1);

    for i = 1:floor(lenM / x)
        Msection = M((i-1)*x+1:i*x,:);
        Mavg = nanmean(Msection);
        M2 = [M2; Mavg];
    end

    % if any bins left...
    if (lenM > (i*x) && ~crop)
        Msection = M(x*i+1:lenM,:);
        if size(Msection,1)>1 
            Mavg = nanmean(Msection);
        else
            Mavg = Msection;
        end
        M2 = [M2; Mavg];
    end
end