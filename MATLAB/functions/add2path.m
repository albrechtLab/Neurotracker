function newpath = add2path(fileLocation, order)
%ADD2PATH   Add a location to MATLAB search path
%   NEWPATH returns a char array of the updated MATLABPATH
%       order: 1-by-1 numeric/double/int/uint value 
%           Determines if location is first or last in the MATLABPATH.
%           (default=1, add new location first)
%           
%   v0.1 V.Kamara   20230124 
   
if ~exist('order', 'var')   order = 1;  end % default: add to beginning of path

oldPath = path;                         % get current path

    if (order == 1)
        newpath = path(fileLocation, oldPath);  % append new location to beginning of matlab path
    else
        newpath = path(oldPath, fileLocation);  % append new location to end of matlab path
    end
    
end
