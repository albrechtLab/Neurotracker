%% Neurotracker Data Import and Analysis

% edit log:
% 6/14/2021 DRA v1.0.0 to produce PDF summary pages for quick evaluation of
%           tracking and experiment performance
% 9/22/2021 DRA v1.1.0 add more explanatory text and realTime function
%           added generic filename parsing, option for no PDF output, and
%           pre-definition of datafolder 
% 6/10/2022 DRA v1.1.1 cleanup
%
% 1/24/2023 VLK v1.1.2 focused on creating an easily functioning version that 
%           can simply be downloaded from the GitHub and ran without much
%           configuration required. Added the path definition immediately
%           beneath (allowing access to PDF generating fncs), removed
%           function definitions from end of file and placed in
%           "./functions", and updated databrowse (formally databrowseS) near 
%           line 174 enforcing a lower bound on indexing values to mitigate 
%           error often generated when creating figure(4) of the summary
%           PDF.
%           

%% Reset matlab env
resetEnvironment = true;
if (resetEnvironment)
    clc; close all; clear all;  % start in clean & clear env
    restoredefaultpath          % test pathdef functionality
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%               DEFINE PATH 2 MATLAB FUNCTIONS FROM GITHUB              %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "neuroTrackerFolder" is the folder holding this script from github.
% "matlabFunctionsFolder" to the functions folder within your local download of the github repository 
% neuroTrackerFolder = 'Z:\Vanessa\Software\GitHub\NeuroTracker\';      % 
% neuroTrackerFolder = pwd;    % repo from github    
% matlabFunctionsFolder = [neuroTrackerFolder,'\functions\'];      % [neuroTrackerFolder,'\functions'];
% requires admin access
% savepath [neuroTrackerFolder,'\functions\']      % save new path in the functions folder

% Add folders to matlab path for proper access
addpath('Z:\Vanessa\Software\GitHub\NeuroTracker;');                    % add repo from github 2 path
addpath('Z:\Vanessa\Software\GitHub\NeuroTracker\functions;');          % add fcns 2 path
addpath('Z:\Vanessa\Software\GitHub\NeuroTracker\functions\PDFtk\bin;') % for pdf funcs
savepath('Z:/Vanessa/Software/GitHub/NeuroTracker/functions/pathdef.m')

% If processing example data, set the stimulus_t to match the experiment parameters.
%   Stimulus was applied for frames 300-350 (showcased by increase in bgnd
%   intensity from fluorescine in solutions).
stimulus_t = [30 35];   % ([startStimulusFrame, stopStimulusFrame]*dt)

%% Define defaults
%   to change, define the variable before running this m-file
%
if ~exist('dt','var') dt = 0.1; end     % frame interval in s (default 0.1s = 10 fps)
if ~exist('baseline_t','var') baseline_t = [0.2 4]; end   % time for calculation of F0 (s)
if ~exist('stimulus_t','var') stimulus_t = [5 15]; end  % time of stimulus on/off (s)
if ~exist('same_animals_per_set','var') same_animals_per_set = true; end  
        % for combined .txt experiment sets set: 
        %   true  if mov/an duplicate represent additional experiments 
        %   false if mov/an duplicate represent additional animals
if ~exist('make_pdf','var') make_pdf = true; end
show_plots = true;
    
fprintf('NOTE: To change settings, cancel (ctrl-C) and change variables in Command Window and rerun.\n');
fprintf('    --e.g., "baseline_t = [0.5 1]" sets baseline F0 time to 0.5-1 s \n');
fprintf('            "elist = [1 5 10 20]" displays experiment trials 1, 5, 10,and 20 in summary plots \n');
fprintf('    --other variables: [dt][baseline_t][stimulus_t][elist][make_pdf][datafolder]\n');
fprintf('\nSettings: \n');        
fprintf('-Time interval (dt): %0.2f s (%0.1f fps) \n',dt, 1/dt);
fprintf('-Baseline F0 (baseline_t):   %0.1f - %0.1f s \n',baseline_t(1),baseline_t(2));
fprintf('-Stimulus time (stimulus_t): %0.1f - %0.1f s \n',stimulus_t(1),stimulus_t(2));
msg = {'false','true'};
fprintf('-Make PDF output (make_pdf): %s \n',msg{(make_pdf>0)+1});
msg = {'animals','experiments'};
fprintf('-Duplicate trials as additional: %s \n\n',msg{(same_animals_per_set>0)+1});

%% Read in neural imaging text files. Ask to reload if data already imported.
%
reloadData = false;
if exist('SmDat','var') 
    ButtonName = questdlg('Reload data?','','Yes', 'No', 'Yes');
    reloadData = strcmp(ButtonName,'Yes');
end
if ~exist('SmDat','var') || reloadData
    if ~exist('datafolder','var') || reloadData
        datafolder = uigetdir(pwd,"Select a folder containing Neurotracker .txt files"); 
    end
    if datafolder==0
        return
    else
        D = dir(fullfile(datafolder,'*.an*.txt'));
        fprintf('Found %d neural trace files. Loading... \n',length(D));
        cd(datafolder);
    end
    
    % Read all NeuroTracker .txt files into SmDat
    clear SmDat
    SmDat(length(D)).name = '';  % preallocate structure
    
    for f = 1:length(D)          % loop though each .txt file
        
        name = fullfile(datafolder,D(f).name);
        SmDat(f).name =   name;
        SmDat(f).data =   textread(SmDat(f).name,'','delimiter',',','headerlines',1);
        SmDat(f).exp =    str2num(name(strfind(name,'mov')+[3:5]));
        SmDat(f).animal = str2num(name(strfind(name,'.an')+[3:4]));
        
        SmDat(f).datenum = datenum(name(strfind(name,'stream_')+[7:22]),'yyyy-mm-dd-HH-MM');
        
        % Parse any other filename labels
        fname = extractAfter(name,'stream_');
        tokfind = bwlabel(isletter(fname));
        for tok = 1:max(tokfind); 
            ix = find(tokfind == tok); iy = find(tokfind == tok+1);
            if ~isempty(iy)  % (if there's a number after a token)
                token{tok} = fname(ix);             % token = list of character labels, if there's a number after it
                value(tok) = str2num(fname((max(ix)+1):(min(iy)-1))); % value is the number after each token
                SmDat(f).(token{tok}) = value(tok);
            end
        end  
        
        if ismember('pattern',name)
            SmDat(f).pattern = str2num(name((strfind(name,'pattern')+3):(strfind(name,'pattern')+6)));
        end
     
        SmDat(f).valve = str2num(name((strfind(name,'valve')+5)));
        if contains(name,'stim')
            % SmDat(f).time = str2num(strrep(name((strfind(string(name),'_mov')+(-5:-1))),'_','.'));
            % SmDat(f).stimulus = str2num(name((strfind(string(name),'stim')+(4:5))));
        else
            SmDat(f).time = [];
            SmDat(f).stimulus = SmDat(f).valve;
        end
        
        if mod(f,round(length(D)/20))==0, fprintf('\b*\n'); end   % display a crude progress bar
        
    end
end
cd(datafolder);

% define filenames and generate or read in a video image frame
disp('Getting video frame...');
[pathname,basename,~] = fileparts(SmDat(1).name); basename = extractBefore(basename,'_mov');
expName = dir('*_settings.txt'); expName = extractBefore(expName.name,'_settings.txt');
videoFile = [extractBefore(SmDat(1).name,'.an'),'.tif'];
imageFile = [extractBefore(SmDat(1).name,'.an'),'.jpg'];
if ~exist(imageFile,'file')
    if ~exist(videoFile,'file')
        [ft,pt] = uigetfile({'*.tif','TIFF stack (*.tif)';'*.jpg','First frame image (*.jpg)'},'Select the first TIFF video stack or JPG image.');
        [~,~,ext] = fileparts(ft);
        if strcmp(ext,'.jpg') 
            imageFile = fullfile(pt,ft);
        elseif strcmp(ext,'.tif') 
            videoFile = fullfile(pt,ft);
        end
    end
    if exist(videoFile,'file')
        img = imread(videoFile,1);
        imwrite(img,imageFile,'BitDepth',12);   
    end
end
img = imread(imageFile);

disp('Organizing data...');

% Parse and align data
raw = struct2mat(3,SmDat,[],{'data'});  % 3-D matrix dimensions: time point x variable x animal
numfiles = size(raw,3);                 % # of .txt files (animals * trials)
numpts = mmax(raw(:,1,:));              % # of recorded timepoints per trial
t = (1:numpts)*dt;                      % vector of real time

% Define time ranges for analysis: baseline (before stim), stimulus, and post-stimulus
%   These are arrays of index values 
baseline_range = find(t > baseline_t(1) & t <= baseline_t(2));
stimulus_range = find(t > stimulus_t(1) & t <= stimulus_t(2));
poststim_range = find(t > stimulus_t(2));

% Handle any overwritten timepoints
All = NaN*ones(numpts,size(raw,2),numfiles);
for ex = 1:numfiles
    valid = find(~isnan(raw(:,1,ex)));
    All(raw(valid,1,ex),:,ex)=raw(valid,:,ex); %this will end up with only the last row if there is more than one row for a frame
end

% Make variables for sorting later: these 1-D vectors are one value per animal/txt file
animal = struct2mat(1,SmDat,[],{'animal'});         % animal number (from filename, an_)
exp = struct2mat(1,SmDat,[],{'exp'});               % experiment/trial number (mov_)
valve = struct2mat(1,SmDat,[],{'valve'});           % if a valve position was defined by filename (e.g. from MVP rotation valve)
expset = cumsum([1;exp(2:end)<exp(1:end-1)]);       % if sets have duplicate mov_ numbers, put them in different sets
% time = struct2mat(1,SmDat,[],{'time'});
% stimulus = struct2mat(1,SmDat,[],{'stimulus'});     % if a stimulus number was defined by filename

% Parse Red Flags in text files (variable column 16)
redFlag = squeeze(All(:,16,:));

% Remove faulty negative values (tracking glitches)
AllSqInt = squeeze(All(:,12,:)); 
AllSqInt(AllSqInt<500)=NaN;  % actually, remove all below 500

% Animal positions
AllXY = squeeze(All(:,2:3,:));                      % get all centroid xy values
XYrange = squeeze(min(AllXY));                      
XYrange = [XYrange; squeeze(max(AllXY)) - XYrange]; % per trial, this is [xmin ymin xrange yrange] in pixels
XY_trial_disp = sqrt(sum(XYrange(3:4,:).^2));         % displacement during trial in pixels


% Renumber animals or expts if multiple individual trackings created duplicates
numsets = max(expset);
oanimal = animal;       % original animal number
oexp = exp;             % original experiment number
if numsets>1
    animalsperset = []; 
    for s=1:numsets 
        animalsperset(s) = max(animal(expset==s))+1;
        exptsperset(s) = max(exp(expset==s));
        if s>1
            animal(expset == s) = animal(expset == s) + sum(animalsperset(1:s-1));
            exp(expset == s) = exp(expset == s) + sum(exptsperset(1:s-1));
        end
    end
end
if (same_animals_per_set) 
    animal = oanimal;       % keep original animal number & adjust expt number
else
    exp = oexp;             % keep original expt number & adjust animal number
end

numanimals = max(animal)+1;
numexp = max(exp);

%% Calculations

disp('Calculating data...');

% Determine range of fluorescence per trial/animal
trial_stats = [nanmin(AllSqInt,[],1); nanmax(AllSqInt,[],1)];
trial_stats = [trial_stats; diff(trial_stats); mean(trial_stats)];

% Calculate baseline fluorescence F0 and normalize F/F0 per trial & animal
baseline = nanmean(AllSqInt(baseline_range,:));
AllSqIntNorm = AllSqInt ./ repmat(baseline,size(AllSqInt,1),1);

% F0 = bcopy;
% F0stimbaseline = baseline;
% % normalized by each individual trial
% AllSqIntNorm = AllSqInt ./ repmat(bcopy,size(AllSqInt,1),1);
% % normalized by last non-voltage trial created above
% AllSqIntF0stimNorm = AllSqInt ./ repmat(baseline,size(AllSqInt,1),1);

% Normalizing Florescence by F0first
F0 = baseline;
F0first = zeros(size(F0));    % F0 of first trial duplicated for all trials
for a = 1:numanimals
    ix = find(animal == a-1);
    F0first(ix) = F0(ix(exp(ix) == min(exp(ix))));
end

% Find raw peaks
extendedStimulusRange = min(stimulus_range):max(stimulus_range)+20;  % add 20 frames to peak search window (stimulus on time extended 20 frames or 2 s)
peakInt = zeros(size(F0));            % peak integrated fluorescence F
peakTime = zeros(size(F0));           % time in seconds at which peak occurs
for i = 1:length(peakInt)
    [peakInt(i), pidx] = nanmax(smooth(AllSqInt(extendedStimulusRange,i)'));
    peakTime(i) = t(extendedStimulusRange(pidx));
end

% % Filter outliers by histogram
% % ADDED line to copy on new matrix
% Flimits = [0.01,0.99]; % exclude lower and upper 1%
% Frange = 0:.1:4;
% h = hist(AllSqIntNorm(:),Frange); ch = cumsum(h) ./ msum(h); 
% hs = hist(AllSqIntNorm(:),Frange); ch = cumsum(hs) ./ msum(hs);% repl AllSqIntNorm
% llim = Frange(find(ch >= Flimits(1),1)-1); ulim = Frange(find(ch >= Flimits(2),1));
% AllSqIntNormFilt = AllSqIntNorm;
% AllSqIntF0stimNormFilt = AllSqIntF0stimNorm;
% AllSqIntNormFilt(AllSqIntNormFilt < llim | AllSqIntNormFilt > ulim) = nan;
% AllSqIntF0stimNormFilt(AllSqIntF0stimNormFilt < llim | AllSqIntF0stimNormFilt > ulim) = nan;
% % did not touch below.
% h = hist(AllSqInt(:),Frange); ch = cumsum(h) ./ msum(h);
% llim = Frange(find(ch >= Flimits(1),1)-1); ulim = Frange(find(ch >= Flimits(2),1));
% AllSqIntFilt = AllSqInt;
% AllSqIntFilt(AllSqIntFilt < llim | AllSqIntFilt > ulim) = NaN;

% Analyze background for stimulus onset (only works if stimulus has fluorescein)
bcg = smoothmat(squeeze(All(:,6,:))); %bcg(1,:)=NaN;% clear 1st frame value
bcg_base = nanmean(bcg(baseline_range,:));          % background intensity before stimulus
bcg_max = max(bcg(stimulus_range,:));               % background intensity during stimulus
bcg_norm = (bcg - repmat(bcg_base,size(bcg,1),1)) ./ repmat(bcg_max - bcg_base,size(bcg,1),1);
        % normalized background: 0 = stim off; 1 = stim on
        
% bcg_animalavg = [];
% for a = 0:(numanimals-1)
%     jj = find(animal == a);
%     bcg_animalavg(:,a+1) = nanmean(bcg_norm(:,jj),2); 
% end
% 
% threshold = 0.2; 
% bcg_start = []; 
% for ii = 1:numanimals
%     bcg_start(ii) = find(bcg_animalavg(1:end-1,ii)<threshold & bcg_animalavg(2:end,ii)>=threshold,1); 
% end
% 
% bcg_start_cycle = (bcg_start - min(bcg_start))*-1;
% bcg_shiftVal = bcg_start_cycle(animal+1);
% 
% bcg_norm_align = [];
% AllSqIntAlign = [];
% for kk = 1:length(animal)
%     bcg_norm_align(:,kk) = circshift(bcg_norm(:,kk),bcg_shiftVal(kk),1);
%     AllSqIntAlign(:,kk) = circshift(AllSqInt(:,kk),bcg_shiftVal(kk),1);
% end



%% -----------SAVE SUMMARY PLOTS-------------
%
if (show_plots) 

    warning('off','MATLAB:linkaxes:RequireDataAxes');
    % (1)
    % Plot selected trials for each animal and label detected peaks and
    % baselines
    %
        fprintf('\n1. Plotting animal summary... ');

        % which trial #s to plot? if many, choose a few (maxcurves)
        % or specify the list via variable 'elist'
        if ~exist('elist','var')
            maxcurves = 4;
            if (numexp <= maxcurves) elist = 1:numexp; 
            else elist = [1 round(logspace(log10(2),log10(numexp),maxcurves-1))]; end
        end
        colors = hsv(max(5,length(elist)))*0.75;

        figure(1); clf;
        spy = ceil(sqrt(numanimals)); spx = ceil(numanimals / spy);
        tl = tiledlayout(spy,spx,'TileSpacing','Compact');
        for a = 1:numanimals
            nexttile;
            for ee = 1:length(elist)
                idx = find(animal == a-1 & exp == elist(ee));
                if ~isempty(idx) 
                    trace = AllSqInt(:,idx)';
                    plot(t,trace,'Color',[0.8 0.8 0.8]); hold on;               % plot raw data
                    plot(t,smooth(trace),'Color',colors(ee,:));                 % plot smoothed data
                    plot(peakTime(idx), peakInt(idx),'o','Color',colors(ee,:)); % circle on peaks
                    plot(0,F0(idx),'*','Color',colors(ee,:));                   % * on baselines
                end
            end
            title(['an ',num2str(a-1)]);
            hold off;
        end

        % make the axis limits nicer...
        try
            nicelims(findobj(gcf,'Type','Axes'),[0 Inf 0 NaN]);
        catch
            linkaxes(get(tl,'Children'),'xy');
        end

        % label axes and print
        hilite(baseline_t,0.05,[0.7 0.7 0.7],'all');       % highlight baseline times
        hilite(stimulus_t,[],[0.8 1 0.9],'all',false);     % highlight stimulus times
        xlabel(tl,'Time(s)'); 
        ylabel(tl,'F (integrated fluorescence)');
        %title(tl,['Peak Int per animal, trials ',mat2str(elist),' ',basename],'Interpreter','none');
        title(tl,['Peak Int per animal, trials ',mat2str(elist),' ',expName],'Interpreter','none');

        if (make_pdf)
            fprintf('saving PDF... ');
            %print([basename,'__PeakInt_perAnimal.pdf'],'-dpdf','-fillpage');
            print([expName,'__PeakInt_perAnimal.pdf'],'-dpdf','-fillpage');
            fprintf('done.');
        end

    % (2)
    % Plot video image and neuron positions

        fprintf('\n2. Plotting animal positions... ');

        figure(2); clf;
        imagesc(imadjust(img));
        cmap = gray;
        colormap(1-cmap(:,1)*[0.5 0.5 0]); 

        hold on;
        for i=1:size(XYrange,2) 
            rectangle('Position',XYrange(:,i),'EdgeColor','r'); 
            % text(XYrange(1,i)+10,XYrange(2,i)+10,num2str(animal(i)),'Color','k','BackgroundColor','w','Margin',1); 
            ix = find(animal == animal(i)); 
            if (i == ix(exp(ix) == min(exp(ix))))  % label just the first trail per animal
                text(XYrange(1,i)+11,XYrange(2,i)+11,num2str(animal(i)),'Color','k');
            end
        end 
        axis ij equal;
        set(gca,'XLim',[0 1024],'YLim',[0 1024]);

        xlabel('X-position (pix)'); 
        ylabel('Y-position (pix)');
        %title(['Animal positions: ',basename],'Interpreter','none');
        title(['Animal positions: ',expName],'Interpreter','none');
        if (make_pdf)
            fprintf('saving PDF... ');
            %print([basename,'__AnimalPos.pdf'],'-dpdf','-fillpage');
            print([expName,'__AnimalPos.pdf'],'-dpdf','-fillpage');
            fprintf('done.');
        end

    %
    % (3)
    % Plot peak, baseline F0 for each animal across trials
    %----------------

        fprintf('\n3. Plotting animal baselines and peaks...');

        figure(3); clf;
        spy = ceil(sqrt(numanimals)); spx = ceil(numanimals / spy);
        tl = tiledlayout(spy,spx,'TileSpacing','Compact');
        for a = 1:numanimals
            nexttile;
            idx = find(animal == a-1); idx = [idx(1);idx];
            plot(exp(idx),[peakInt(idx);F0(idx);F0first(idx)]','.-');
            title(['an ',num2str(a-1)]);
            if (a==1) legend({'peak','F_0','F_0(first)'},'Location','none','Box','off'); end
        end
        % make the axis limits nice...
        try
            nicelims(findobj(gcf,'Type','Axes'),[0 numexpts+1 0 NaN]);
        catch
            linkaxes(get(tl,'Children'),'xy');
        end
        xlabel(tl,'Trial #');
        ylabel(tl,'F (integrated fluorescence)');
        %title(tl,['Trial Peak, F0, and F0(first) per animal: ',basename],'Interpreter','none');
        title(tl,['Trial Peak, F0, and F0(first) per animal: ',expName],'Interpreter','none');
        if (make_pdf)
            fprintf('saving PDF... ');
            %print([basename,'__Peak_F0_perAnimalperTrial.pdf'],'-dpdf','-fillpage');
            print([expName,'__Peak_F0_perAnimalperTrial.pdf'],'-dpdf','-fillpage');
            fprintf('done.');
        end

%     %
%     % (4)
%     % Plot background (stimulus check)
% 
        fprintf('\n4. Plotting background (stimulus check)...');

        figure(4); clf; tl = tiledlayout(2,1);
        ix = find(t>0.5); % exclude first 0.5s for proper scaling
        h = databrowse(bcg(ix,:)',t(ix),exp);       
        figure(4); nexttile; 
        copyobj(get(h(1),'Children'),gca);
        title('Stimulus background by experiment/trial');

        % plot normalized background fluorescence per animal
        h = databrowse(bcg(ix,:)',t(ix),animal);
        figure(4); nexttile;
        copyobj(get(h(1),'Children'),gca);
        title('Stimulus background by animal');
        hilite(baseline_t,0.05,[0.7 0.7 0.7],'all');       % highlight baseline times
        hilite(stimulus_t,[],[0.8 1 0.9],'all',false);     % highlight stimulus times

        xlabel(tl,'Time(s)'); 
        ylabel(tl,'Background');
        %title(tl,['Stimulus check: ',basename],'Interpreter','none');
        title(tl,['Stimulus check: ',expName],'Interpreter','none');
        if (make_pdf)
            fprintf('saving PDF... ');
            %print([basename,'__StimulusCheck.pdf'],'-dpdf','-fillpage');
            print([expName,'__StimulusCheck.pdf'],'-dpdf','-fillpage');
            fprintf('done.');
        end

        disp(' All done.');

    % summarize PDFs into a single file
    % (this function needs pdftk.exe found via pathdef at beginning
    if (make_pdf)
        %onepdf(dir('*__*.pdf'),[basename,'_Summary.pdf'],true);
        onepdf(dir('*__*.pdf'),[expName,'_Summary.pdf'],true);
    end
% 
end    
    
%% Display data browser
% databrowse(AllSqIntNorm',t,exp);
% 
% disp('View data using these examples:');
% disp('databrowse(AllSqIntNorm'',t,exp);    % group by experiment');
% disp('databrowse(AllSqIntNorm'',t,animal); % group by animal');
% disp('idx = find(ismember(animal,[1,5])); databrowse(AllSqIntNorm(:,idx)'',t,expset(idx)); % group by experiment set, only animals 1 and 5');
% 
% RT = realTimeData(AllSqIntNorm,exp,animal);
%databrowse(RT.mat',[],RT.animal);


