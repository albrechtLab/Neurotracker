function output = anova1multicompare(datavector, groupvector, plist, comparisontype)

%-----------------------------------------------------------------------
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