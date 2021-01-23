%% Neurotracker Data Import and Analysis

% Read in neural imaging text files. Ask to reload if data already imported
if exist('SmDat','var') 
    ButtonName = questdlg('Reload data?','','Yes', 'No', 'Yes');
end
if ~exist('SmDat','var') || strcmp(ButtonName,'Yes')
    p = uigetdir(pwd,"Select a folder containing Neurotracker .txt files");
    if p==0
        return
    else
        D = dir(fullfile(p,'*.an*.txt'));
        fprintf('Found %d neural trace files.  \n',length(D));
        cd(p);
    end
    
    % Read all NeuroTracker .txt files into SmDat
    clear SmDat
    SmDat(length(D)).name = '';  % preallocate structure
    
    for f = 1:length(D) 
        name = fullfile(p,D(f).name);
        SmDat(f).name = name;
        SmDat(f).data = textread(SmDat(f).name,'','delimiter',',','headerlines',1);
        SmDat(f).exp = str2num(name((strfind(name,'mov')+3):(strfind(name,'mov')+5)));
        %SmDat(f).pattern = str2num(name((strfind(name,'pattern')+3):(strfind(name,'pattern')+6)));
        if ismember('pattern',name)
            SmDat(f).pattern = str2num(name((strfind(name,'pattern')+3):(strfind(name,'pattern')+6)));
        end
        
        SmDat(f).animal = str2num(name((strfind(name,'.an')+3):(strfind(name,'.an')+4)));
        SmDat(f).valve = str2num(name((strfind(name,'valve')+5)));
        %SmDat(f).pattern = str2num(name((strfind(name,'pattern')+7)));
        if contains(name,'stim')
            SmDat(f).time = str2num(strrepl(name((strfind(name,'_mov')+(-5:-1))),'_','.'));
        else
            SmDat(f).time = [];
        end
        if contains(name,'stim')
            SmDat(f).stimulus = str2num(name((strfind(name,'stim')+(4:5))));
        else
            SmDat(f).stimulus = SmDat(f).valve;
        end
        if mod(f,round(length(D)/20))==0, fprintf('\b*\n'); end   % display a crude progress bar
    end
end

% Parse and align data
raw = struct2mat(3,SmDat,[],{'data'});   % 3-D matrix dimensions: time point x variable x animal
numfiles = size(raw,3);
numpts = mmax(raw(:,1,:));

% Handle any overwritten timepoints
All = NaN*ones(numpts,size(raw,2),numfiles);
for ex = 1:numfiles
    valid = find(~isnan(raw(:,1,ex)));
    All(raw(valid,1,ex),:,ex)=raw(valid,:,ex); %this will end up with only the last row if there is more than one row for a frame
end

% Make variables for sorting later: these 1-D vectors are one value per animal/txt file
animal = struct2mat(1,SmDat,[],{'animal'});
exp = struct2mat(1,SmDat,[],{'exp'});
%pattern = struct2mat(1,SmDat,[],{'pattern'});
%exp = All(2,15,:);
valve = struct2mat(1,SmDat,[],{'valve'});
%pattern = struct2mat(1,SmDat,[],{'pattern'});
expset = cumsum([1;exp(2:end)<exp(1:end-1)]);
% time = struct2mat(1,SmDat,[],{'time'});
stimulus = struct2mat(1,SmDat,[],{'stimulus'});

% Parse Red Flags in text files (variable column 16)
redFlag = squeeze(All(:,16,:));

% Remove faulty negative values (tracking glitches)
AllSqInt = squeeze(All(:,12,:)); 
AllSqInt(AllSqInt<500)=NaN;  % actually, remove all below 500

% Renumber animals if multiple individual trackings created duplicates
numsets = max(expset);
if numsets>1
    oanimal = animal;
    animalsperset = []; 
    for s=1:numsets 
        animalsperset(s) = max(animal(expset==s))+1;
        if s>1
            animal(expset == s) = animal(expset == s) + sum(animalsperset(1:s-1));
        end
    end
end
numanimals = max(animal)+1;

% Calculate baseline fluorescence F0 and normalize F/F0
baseline = nanmean(AllSqInt(1:40,:));
AllSqIntNorm = AllSqInt ./ repmat(baseline,size(AllSqInt,1),1);

% Filter outliers by histogram
Flimits = [0.01,0.99]; % exclude lower and upper 1%
Frange = 0:.1:4;
h = hist(AllSqIntNorm(:),Frange); ch = cumsum(h) ./ msum(h);
llim = Frange(find(ch >= Flimits(1),1)-1); ulim = Frange(find(ch >= Flimits(2),1));
AllSqIntNormFilt = AllSqIntNorm;
AllSqIntNormFilt(AllSqIntNormFilt < llim | AllSqIntNormFilt > ulim) = NaN;

h = hist(AllSqInt(:),Frange); ch = cumsum(h) ./ msum(h);
llim = Frange(find(ch >= Flimits(1),1)-1); ulim = Frange(find(ch >= Flimits(2),1));
AllSqIntFilt = AllSqInt;
AllSqIntFilt(AllSqIntFilt < llim | AllSqIntFilt > ulim) = NaN;

% Display data browser
databrowseS(AllSqIntNorm');

% figure(1); 
% for a = 0:max(animal) 
%     subplotarray(numanimals,a+1); 
%     imagesc(AllSqIntNorm(:,animal == a)'); 
%     set(gca,'CLim',[.5 4]); 
%     title(a); 
% end

figure(3); clf; 
for a = 0:max(animal) 
    lineartrace = reshape(AllSqIntNorm(:,animal == a),[],1);
    smoothtrace = smooth(lineartrace,10);
    smoothtrace(isnan(lineartrace)) = NaN;
    hold on; 
    plot(smoothtrace + (max(animal)-1-a)*4); 
end; hold off


bcg = smoothmat(squeeze(All(:,6,:)));
bcg(1,:)=NaN; % clear 1st frame value
baseline_bcg = nanmean(bcg(5:40,:));
baseline_init = 0;
bcg_max = max(bcg(50:120,:));
bcg_norm2 = (bcg - repmat(baseline_bcg,size(bcg,1),1)) ./ repmat(bcg_max - baseline_bcg,size(bcg,1),1);

animal2 = animal+1;
anrange = 1:max(animal2);


bcg_animalaverage = [];
for a = anrange
    jj = find(animal2 == a);
    bcg_animalaverage(:,a) = nanmean(bcg_norm2(:,jj),2); 
end

threshold = 0.2; 
bcg_start = []; 
for ii = anrange
    bcg_start(ii) = find(bcg_animalaverage(1:end-1,ii)<threshold & bcg_animalaverage(2:end,ii)>=threshold,1); 
end

bcg_start_cycle = (bcg_start - min(bcg_start))*-1;
bcg_shiftVal = bcg_start_cycle(animal2);

bcg_norm_align = [];
AllSqIntNormAlign = [];
AllSqIntNormFiltAlign = [];
for kk = 1:length(animal);
    bcg_norm_align(:,kk) = circshift(bcg_norm2(:,kk),bcg_shiftVal(kk),1);
    AllSqIntNormAlign(:,kk) = circshift(AllSqIntNorm(:,kk),bcg_shiftVal(kk),1);
    AllSqIntNormFiltAlign(:,kk) = circshift(AllSqIntNormFilt(:,kk),bcg_shiftVal(kk),1);
end

ASINAmax = max(AllSqIntNormAlign);
ASINAmr = []; %ASINA max reorder
for a = 0:max(animal) 
    tracemax = reshape(ASINAmax(:,animal == a),[],1);
    ASINAmr = [ASINAmr, tracemax];
end


lineartracevals = [];
for a = 0:max(animal) 
    lineartracevals = [lineartracevals, reshape(AllSqIntNorm(:,animal == a),[],1);];
end

lineartraceRedFlag = [];
for a = 0:max(animal) 
    lineartraceRedFlag = [lineartraceRedFlag, reshape(redFlag(:,animal == a),[],1);];
end


%Remove red flagged animal traces from lineartracevals dataset
redFlagIdx = find(any(lineartraceRedFlag == 1));

ASINAmr(:,redFlagIdx) = [];
lineartracevals(:,redFlagIdx) = [];

clear redFlag lineartraceRedFlag redFlagIdx

%databrowse(ASINAmr')



%% HELPER FUNCTIONS

%% struct2mat 
%   pulls structure information into matrices

function output = struct2mat(dim,structure,index,fields)

    if ~iscell(fields) fields = {fields}; end
    if isempty(index) index = 1:numel(structure); end

    output = [];
    for i = 1:length(index)
        substructure = structure(index(i));
        for j = 1:length(fields)
            substructure = getfield(substructure,char(fields(j)));
        end

        output = safecat(dim,output,substructure);
    end
end
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
