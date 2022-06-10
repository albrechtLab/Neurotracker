
function handle = databrowseS(mat,x,ygroup,defaultsettings)
%
% function handle = databrowseS(mat,x,ygroup,defaultsettings)
% 
%          mat: m-by-n data matrix
%            x: 1-by-n optional labels for x-axis
%       ygroup: m-by-1 list of group nnames for each row in mat
%                 (can be numeric or cell)
%
%    defaultsettings: cell array containing values of all user objects
%           e.g. {1:length(groups),0,1,0,1,1,1,0}
%           (1) groups to display, (2) individual traces, (3) mean,
%           (4) SD, (5) SEM, (6) shaded error, (7) x-binning, 
%           (8) show stats, (9) sort group order, (10) sort matrix rows     
%
%       Note: can transfer figure elsewhere by: 
%             h = databrowse(...); copyobj(get(h(1),'Children'),gca);
%
%   v1.1 D.Albrecht 20100427
%   v1.5 D.Albrecht 20210111


% find or make DataBrowse figure
figH = findobj(0,'tag','DataBrowse');
if isempty(figH)
    figH = figure(25); set(figH,'Tag','DataBrowse','Name','DataBrowse','NumberTitle','off');
end
figure(figH);

% find existing data if none given
if ~exist('mat')
    Data = get(figH,'UserData');
    if ~isempty(Data)
        mat = Data.mat;
        x = Data.x;
        ygroup = Data.ygroup;
        h = Data.handles;
        sig = Data.sig; 
        output = Data.output;
    end
end
if ~exist('x') || isempty(x) x = 1:size(mat,2); end
if ~exist('ygroup') || isempty(ygroup) ygroup = ones(size(mat,1),1); end
if ~exist('sig') sig = []; end
if ~exist('output') output = []; end

ds = {1:length(unique(ygroup)),0,1,0,1,1,1,0,1,0};
if ~exist('defaultsettings') 
    defaultsettings = ds;
else
    deflen = length(defaultsettings);
    if deflen < length(ds)
        defaultsettings(deflen+1:length(ds)) = ds(deflen+1:length(ds));
    end
end

% default settings
if ~exist('Data','var') || ~isfield(Data,'settings')
    Data.settings = defaultsettings;
else
    Data.settings = get(h(3:12),'Value');
end
values = [1:6,8,10,15,20,25,30,-1,-2];

sorted = false;  % check if either sorting method is used. this is needed as stats are calcuated with a sorted group order

% sort matrix rows
if cell2mat(Data.settings(10))
    [~,idx]=sort(ygroup); 
    sorted = true;
else
    idx = 1:length(ygroup);
end

% get group info
[a,b,c] = unique(ygroup(idx)); [~,e] = sort(b); [~,grord] = sort(e);
if cell2mat(Data.settings(9))       % sort groups
    groups = a; groupindex = c; 
    sorted = true;            
else                                % don't sort groups
    groups = a(e); groupindex = grord(c); 
end
if sorted statord = 1:max(e); else statord = e; end % needed for stats, which calcuates using sorted groups!

% set up figure if needed
if ~exist('h')
    h(1) = subplot('Position',[.1 .55 .7 .4]);
    h(2) = subplot('Position',[.1 .1 .7 .35]);

    h(3) = uicontrol('Style', 'listbox',...
           'Tag','grouplist',...
           'Units','normalized',...
           'Position', [0.85, 0.55, 0.1, 0.4],...
           'String', groups,...
           'Value', cell2mat(Data.settings(1)),...
           'Min',1,'Max',length(groups)+1,...
           'Tooltip','Select data groups to graph. Use <shift> or <ctrl> to multiselect.',...
           'CallBack', 'databrowseS;');  
       
    h(4) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.4, 0.1, 0.05],...
           'String', 'Individual',...
           'Value', cell2mat(Data.settings(2)),...
           'Tooltip','Show or hide individual data',...
           'CallBack', 'databrowseS;');
       
    h(5) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.35, 0.1, 0.05],...
           'String', 'Mean',...
           'Value', cell2mat(Data.settings(3)),...
           'Tooltip','Show or hide group mean data',...
           'CallBack', 'databrowseS;');
       
    h(6) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.3, 0.1, 0.05],...
           'String', 'SD',...
           'Value', cell2mat(Data.settings(4)),...
           'Tooltip','Show or hide group standard deviation (SD)',...
           'CallBack', 'databrowseS;');
       
    h(7) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.25, 0.1, 0.05],...
           'String', 'SEM',...
           'Value', cell2mat(Data.settings(5)),...
           'Tooltip','Show or hide group standard eror of the mean (SEM)',...
           'CallBack', 'databrowseS;');
       
    h(8) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.2, 0.1, 0.05],...
           'String', 'Fill',...
           'Value', cell2mat(Data.settings(6)),...
           'Tooltip','Show SD/SEM as shading (or error bars)',...
           'CallBack', 'databrowseS;');
   
    h(9) = uicontrol('Style', 'popupmenu',...
           'Units','normalized',...
           'String',sprintf('%d|',values),...
           'Position', [0.85, 0.075, 0.1, 0.05],...
           'Value', cell2mat(Data.settings(7)),...
           'Tooltip','Bin data over time by the selected factor. Negative values remove outliers.',...
           'CallBack', 'databrowseS;');

    h(10) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.15, 0.1, 0.05],...
           'String','Stats',...           
           'Value', cell2mat(Data.settings(8)),...
           'Tooltip','Show or hide statistics\n*CAREFUL* Sorting may introduce errors, be sure to double check',...
           'CallBack', 'databrowseS;');
       
    h(11) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.85, 0.95, 0.1, 0.025],...
           'String','Sort groups',...           
           'Value', cell2mat(Data.settings(9)),...
           'Tooltip','Sort data group list',...
           'CallBack', 'databrowseS;');
       
    h(12) = uicontrol('Style', 'togglebutton',...
           'Units','normalized',...
           'Position', [0.05, 0.95, 0.05, 0.025],...
           'String','Sort',...           
           'Value', cell2mat(Data.settings(10)),...
           'Tooltip','Sort matrix data rows by group',...
           'CallBack', 'databrowseS;');
       
    % set colors       
    ng = max(4,length(groups)); sep = floor(254/(ng-1));
    a = jet(sep*(ng) + 1); a = a(1:sep:255,:);
    set(h(2),'ColorOrder',a);
else
    set(h(3),'String',groups);
end

% X-binning
frameavg = values(cell2mat(Data.settings(7)));
if frameavg > 1
    mat2 = binaverage(mat(idx,:)',frameavg,0)';
    x2 = x(1:frameavg:end);
else
    mat2 = mat(idx,:);
    x2 = x;
end

% plot heatmap
xn = x2; if iscell(x2) xn = 1:length(x); end
subplot(h(1)); 
imagesc(xn,1:size(mat2,1),mat2); xlabel('Time');
set(gca,'YTick',1:size(mat2,1),'YTickLabel',ygroup(idx));
colormap([.7 .7 .7; jet(254)]);
hc = colorbar;
cbpos = get(hc,'Position'); cbpos(1)=0.81; cbpos(3)=0.01; set(hc,'Position',cbpos);
if iscell(x) set(gca,'XTick',xn,'XTickLabel',num2str(cell2mat(x))); end
if cell2mat(Data.settings(10))
    ticks = ygroup(idx);
    if ~iscell(ticks) ticks = cellstr(num2str(ticks(:))); end
    ticks(find(strcmp(ticks(2:end),ticks(1:end-1)))+1)={'`'};   % simplify matrix row labels
    set(gca,'YTickLabel',ticks);
end

% get ui settings
groupsel = cell2mat(Data.settings(1));
showind = cell2mat(Data.settings(2));
showmean = cell2mat(Data.settings(3));
showsd = cell2mat(Data.settings(4));
showsem = cell2mat(Data.settings(5));
plotfilled = cell2mat(Data.settings(6));
showstats = cell2mat(Data.settings(8));

% Calculate statistics
if showstats && (isempty(sig) || size(mat2,2) ~= size(sig,2))
    %allexps = iselement(groupindex,groupsel);
    tic
    clear sig;
    for i = 1:size(mat2,2) 
        %output = anova1multicompare(mat2(:,i),ygroup(idx));
        output = anova1multicompare(mat2(:,i),groupindex);
        %output = anova1multicompare(mat2(:,i),ygroup);
        sig(:,i) = output.stats(:,4); 
    end
    %sig = sig(grord,:);
    toc
end

% plot traces
subplot(h(2)); cla; hold on;

ydat = [];
for g = 1:length(groupsel)
    exps = find(groupindex == groupsel(g));

    % TEST -- remove outliers
    sigma = 3;
    if frameavg < 1
        
        % get uncorrected mean/sd
        if length(exps)>1
            ymean = nanmean(mat(exps,:));
            ysd = nanstd(mat(exps,:));
        else
            ymean = mat(exps,:);
            ysd = zeros(size(ymean));
        end

        temp = mat(exps,:);
        outlier = (temp > repmat(ymean+sigma*ysd,size(temp,1),1)) | (temp < repmat(ymean-sigma*ysd,size(temp,1),1));
        temp(find(outlier))=NaN;
        mat3(exps,:) = temp;
        
        if frameavg < -1
            mat2 = binaverage(mat3',-frameavg,0)';
            x2 = x(1:-frameavg:end);
        else
            mat2 = mat3;
            x2 = x;
        end

        xn = x2; if iscell(x2) xn = 1:length(x); end
    end

    if length(exps)>1
        ymean = nanmean(mat2(exps,:));
        ysd = nanstd(mat2(exps,:));
        ysem = ysd / sqrt(length(exps));
    else
        ymean = mat2(exps,:);
        ysd = zeros(size(ymean)); ysem = ysd;
    end

    colors = get(gca,'ColorOrder');
    color = colors(rem(groupsel(g)-1,size(colors,1))+1,:);
    
    if showind
        plot(xn,mat2(exps,:),'Color',.6+.2*color);
    end
    
    if showsd
        if plotfilled
            jbfill(xn, ymean-ysd, ymean+ysd, .8+.2*color, .8+.2*color, 1, 1);
        else
            errorbar(xn, ymean, ysd, 'Color', .6+.4*color);
        end
    end
    hold on;
    
    if showsem
        if plotfilled
            jbfill(xn, ymean-ysem, ymean+ysem, .6+.4*color, .6+.4*color, 1, 1);
        else
            errorbar(xn, ymean, ysem, 'Color', .8*color);
        end
    end
    hold on;

    if showmean
        plot(xn,ymean,'LineWidth', 2,'Color', .7*color, 'Tag','mean');
    end
    
    if showstats & length(groupsel)>1 & g>1
        pv = find(output.stats(:,1) == groupsel(1) & iselement(output.stats(:,2),statord(groupsel(g))));
        scatter(xn,ymean,20*(sig(pv,:)+0.1).^2,.5*color)
    end
    
    ydat(g,:,1) = ymean;
    ydat(g,:,2) = ysd;
    ydat(g,:,3) = ysem;
end
if iscell(x) set(gca,'XTick',xn,'XTickLabel',num2str(cell2mat(x))); end
xlabel('Time');
grid on;

%set(findobj(gca,'type','patch'),'EdgeAlpha',0,'FaceAlpha',0.6);

% set shading to background
ch = get(gca,'Children');
ch2 = [findobj(ch,'Type','line'); findobj(ch,'Type','scatter'); findobj(ch,'Type','patch'); findobj(ch,'Type','errorbar')];
%set(gca,'Children',sort(ch,'descend'));
set(gca,'Children',ch2);

% Add legend
legh = findobj(gca,'Tag','mean');
if ~isempty(legh)
    if isnumeric(groups)
        if diff(size(groups))>0 groups = groups'; end
        grouplabel = cellstr(num2str(groups)); 
    else
        grouplabel = groups;
    end
    if true  % add (n= ...) to legend
        for g = 1:length(groups)
            n = length(find(groupindex == g));
            grouplabel{g} = [grouplabel{g},' (n=',num2str(n),')'];
        end
    end
    legend(flipud(legh), grouplabel(groupsel));
end

% deposit data in figure
Data.mat = mat;
Data.x = x;
Data.ygroup = ygroup;
Data.handles = h;
Data.sig = sig;
Data.output = output;
Data.ydata = ydat;
set(figH,'UserData',Data);


handle = h([2,1]);

end

%-----------------------------------------------------------------------
function output = anova1multicompare(datavector, groupvector, plist, comparisontype)
% USAGE: output = anova1multicompare(datavector, groupvector, plist, comparisontype)
%
%       plist is list of p-values to test (default is : [0.05 0.01 0.001 0.0001]
%       comparisontype is 'tukey-kramer' (default), 'bonferroni', etc.
%

if nargin < 4 comparisontype = 'tukey-kramer'; end
if nargin < 3 plist = [0.05 0.01 0.001 0.0001]; end

plist = [1,sort(plist,'descend')];

ds = size(datavector);
gs = size(groupvector);
if ~all(ds == gs)
    if all(ds == fliplr(gs))
        gs = gs';
    else
        error('datavector and groupvector must have same size');
        return
    end
end
    
[p,t,st] = anova1(datavector,groupvector,'off');

output.anovap = p;
output.anovat = t;
output.anovast = st;

n = 0; for i = 1:length(st.means)-1; n = n+i; end
pvali = ones(n,1);
for i = 2:length(plist)

    [c,m,h,gnames] = multcompare(st,'alpha',plist(i),'ctype',comparisontype,'display','off');

    ptest = ~xor(c(:,3)>0, c(:,5)>0);  % true if comparison p-value at least plist(i)

    pvali = pvali + ptest;

    output.multcomp(i-1).alpha = plist(i);
    output.multcomp(i-1).c = c;
end

pval = plist(pvali);
ns = find(pvali == 1);

output.stats = [c(:,[1 2 4]), pvali-1, pval'];
output.statcell = [gnames(c(:,1:2)), num2cell(c(:,4)), num2cell(pval')];
output.statcell(ns,4) = {'n.s.'};
output.groups = gnames;

end

%% jbfill
%   plot filled graph polygons

function[fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency,stair)
%USAGE: [fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency,stair)
%This function will fill a region with a color between the two vectors provided
%using the Matlab fill command.
%
%fillhandle is the returned handle to the filled region in the plot.
%xpoints= The horizontal data points (ie frequencies). Note length(Upper)
%         must equal Length(lower)and must equal length(xpoints)!
%upper = the upper curve values (data can be less than lower)
%lower = the lower curve values (data can be more than upper)
%color = the color of the filled area 
%edge  = the color around the edge of the filled area
%add   = a flag to add to the current plot or make a new one.
%transparency is a value ranging from 1 for opaque to 0 for invisible for
%the filled color only.
%
%John A. Bockstege November 2006;
%Example:
%     a=rand(1,20);%Vector of random data
%     b=a+2*rand(1,20);%2nd vector of data points;
%     x=1:20;%horizontal vector
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1))
%     grid on
%     legend('Datr')
    if nargin<8;stair=0;end %default is to have smooth curve, not stairs
    if nargin<7;transparency=.5;end %default is to have a transparency of .5
    if nargin<6;add=1;end     %default is to add to current plot
    if nargin<5;edge='k';end  %dfault edge color is black
    if nargin<4;color='b';end %default color is blue

    %-----------------
    % DA added
    upper(isnan(upper))=0;
    lower(isnan(lower))=0;
    if edge == -1
        edge = 'k';
        noedge = 1;
    else
        noedge = 0;
    end
    if stair
        if diff(size(xpoints))<1 xpoints = xpoints'; end
        if diff(size(upper))<1 upper = upper'; end
        if diff(size(lower))<1 lower = lower'; end
        dx = mean(diff(xpoints));
        xpoints = reshape(addortho([-dx/2; dx/2],xpoints),[],1)';
        upper = reshape(repmat(upper,2,1),[],1)';
        lower = reshape(repmat(lower,2,1),[],1)';
    end
    %------------------

    if length(upper)==length(lower) && length(lower)==length(xpoints)
        msg='';
        filled=[upper,fliplr(lower)];
        xpoints=[xpoints,fliplr(xpoints)];
        if add
            hold on
        end
        fillhandle=fill(xpoints,filled,color);%plot the data
        set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color

        %-----------
        % DA added
        if noedge
            set(fillhandle,'LineStyle','none');
        end

        if add
            hold off
        end
        set(gca,'Layer','top');
    else
        msg='Error: Must use the same number of points in each vector';
    end
end

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


function result = iselement(M, condition)

    valid = zeros(size(M));
    for c = 1:length(condition)
        valid = valid | (M == condition(c));
    end

    result = valid;
end


function hilite(x0,y0,col,axlist,replace)
% function hilite(x,y,col,axlist,replace)
%       x is Nx2 matrix of [start end] x-pairs for rectangular highlight
%       y is Nx2 matrix of [start end] y-pairs for rectangular highlight
%           either x or y can be [] to span axis limits
%       color is 1x3 RGB color value (or [] = default)
%       axlist is list of axes ([] is current axis, or 'all' for all axes in current figure)
%       replace if true, clears prior highlighting (default = true)
%
% v1.1 2010-04-08
% Dirk Albrecht 

    if ~exist('replace','var') || isempty(replace) replace = true; end
    if ~exist('axlist','var') || isempty(axlist) axlist = gca; end
    if ischar(axlist) axlist = findobj(gcf,'type','axes'); end

    if ~exist('col','var') || isempty(col) col = [0.8 0.8 0.8]; end

    if ~exist('x0','var') x0 = []; end 
    if ~exist('y0','var') y0 = []; end 

    if ~isempty(isnan(x0)) && any(isnan(x0)) axlist = []; end  % no hilighting if NaN
    if ~isempty(isnan(y0)) && any(isnan(y0)) axlist = []; end  % no hilighting if NaN

    for a = 1:length(axlist)
        ax = axlist(a);
        axlim = [get(ax,'XLim'), get(ax,'YLim')];

        if isempty(x0) x = axlim(1:2); else x = x0; end
        if isempty(y0) y = axlim(3:4); else y = y0; end

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
    end
end