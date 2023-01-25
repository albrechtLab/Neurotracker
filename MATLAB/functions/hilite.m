%% hilite
%

% Highlights defined rectangular region(s) in axis window(s).
%
% USAGE:
%    hilite(x,y,col,axlist,replace)
%
%       x is Nx2 matrix of [start end] x-pairs for rectangular highlight
%       y is Nx2 matrix of [start end] y-pairs for rectangular highlight
%           either x or y can be [] to span axis limits
%           or a number 0-1.0 representing thefraction to shade from left
%           or from top
%       color is 1x3 RGB color value (or [] = default)
%       axlist is list of axes ([] is current axis, or 'all' for all axes in current figure)
%       replace if true, clears prior highlighting (default = true)

%---------------------------- 
% Dirk Albrecht 
% V1.3 27-May-2021 added partial shading option
%---------------------------- 

function hilite(x0,y0,col,axlist,replace)

    if ~exist('replace','var') || isempty(replace) replace = true; end
    if ~exist('axlist','var') || isempty(axlist) axlist = gca; end
    if ischar(axlist) axlist = findobj(gcf,'type','axes'); end
        
    if ~exist('col','var') || isempty(col) col = [0.8 0.8 0.8]; end
    
    if ~exist('x0','var') x0 = []; end 
    if ~exist('y0','var') y0 = []; end 
    
    for a = 1:length(axlist)
        ax = axlist(a);
        axlim = [get(ax,'XLim'), get(ax,'YLim')];
        
        if isempty(x0) x = axlim(1:2); 
            elseif (length(x0)==1) x = [axlim(1) axlim(1)+x0*diff(axlim(1:2))];
            else x = x0; end
        if isempty(y0) y = axlim(3:4); 
            elseif (length(y0)==1) y = [axlim(4)-y0*diff(axlim(3:4)) axlim(4)];
            else y = y0; end
        if isnan(x0) x = []; end
            
        % clear prior highlighing
        if replace
            h = findobj(ax,'Tag','hilite');
            if ~isempty(h) delete(h); end
        end
        
        axchildren = get(ax,'Children');
    
        h = [];
        for i = 1:size(x,1)
            for j = 1:size(y,1)
                h(i,j) = patch(x(i,[1 2 2 1 1]),y(j,[1 1 2 2 1]),col);
            end
        end
    
        set(h,'LineStyle','none','Parent',ax,'Tag','hilite');
    
        % order highlighting to back
        newhandles = reshape(h,[],1);
        set(ax,'Children',[axchildren; newhandles]);
        set(ax,'Layer','top');
    end
      
end