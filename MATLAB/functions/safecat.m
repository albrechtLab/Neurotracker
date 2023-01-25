%% safecat 
%   allows for error-free concatenation of matrices of different size

function output = safecat(dim,A,B,blank)
% output = safecat(dim,A,B,blank)
%           like cat function but works with matrices of different size.  
%           Now N-dim arrays supported also.
%           blank is default [NaN]... undefined elements set to blank

    if nargin < 4 blank = NaN; end

    maxdims = max([ndims(A), ndims(B), dim]);
    sizes = [size(A), repmat(1,1,maxdims - ndims(A)); ...
             size(B), repmat(1,1,maxdims - ndims(B))];

    % test for empty matrix
    if isempty(A) 
        output = B;
    elseif isempty(B)
        output = A;
    else
        maxsize = max(sizes);

        if any(sizes(1,:) < maxsize)
            dimexpand = find(sizes(1,:) < maxsize);
            dimexpand = dimexpand(find(dimexpand ~= dim)); % don't expand on 'dim' dimension
            for i = 1:length(dimexpand)
                S.type = '()'; S.subs = {};
                for j = 1:length(maxsize); S.subs{j}=1:maxsize(j); end
                S.subs{dim} = 1:sizes(1,dim);
                S.subs{dimexpand(i)} = (sizes(1,dimexpand(i))+1):maxsize(dimexpand(i));
                A = subsasgn(A,S,blank); 
            end
        end

        if any(sizes(2,:) < maxsize)
            dimexpand = find(sizes(2,:) < maxsize);
            dimexpand = dimexpand(find(dimexpand ~= dim)); % don't expand on 'dim' dimension
            for i = 1:length(dimexpand)
                S.type = '()'; S.subs = {};
                for j = 1:length(maxsize); S.subs{j}=1:maxsize(j); end
                S.subs{dim} = 1:sizes(2,dim);
                S.subs{dimexpand(i)} = (sizes(2,dimexpand(i))+1):maxsize(dimexpand(i));
                B = subsasgn(B,S,blank); 
            end
        end

        output = cat(dim,A,B);
    end
end