function l0_epoch_feedback(subject_file, pathStruct, binlist)

%% Epoch data and detect artifact in epochs using EEGLAB and ERPLAB artifact detection procedures

% Last modified by Hause Lin 10-05-19 1:31 PM hauselin@gmail.com

%% Set up paths and parameters

% subject = '1'; % subject id for debugging

%% get the subject file name without the extension
[~, filetoread_noExt, ~] = fileparts(subject_file);

%% define the time format to use when saving error/log files
tformat = 'dd-mmm-yyyy - HH MM SS AM';

% paths

PATHS = struct();
PATHS.cwd = pwd;
PATHS.parentdir = PATHS.cwd(1:max(strfind(PATHS.cwd,filesep)));
PATHS.dataDir = [strip(pathStruct.matdatadir, 'right', filesep) 'Preprocessed_ICA_cleaned']; % directory containing ICA cleaned subject data
PATHS.subjdatafile = [filetoread_noExt '_continuous_icacleaned.mat'];
PATHS.epochdir = fullfile(PATHS.cwd,'Epochs');
PATHS.binlister = binlist;
PATHS.log = fullfile(PATHS.cwd,'Log');
PATHS.artifact = fullfile(PATHS.cwd,'Artifacts');
PATHS.designmat = fullfile(PATHS.cwd,'DesignMatrix');

if ~exist(PATHS.binlister, 'file')
    disp('The specified binlist file path does not exist.')
    return
end

% epoching configuration parameters
CFG = struct();
CFG.event0 = 'stimulus'; % timelocking event
CFG.epochRange = [-2500 2500]; % epoch duration around timelocking event
CFG.epochBaselineCorrect = [-200 0]; % baseline correction window ms
CFG.rejecttrialpercent = 25; % mark subjects with more than this % of trials with artifacts
CFG.epochSave = true; % save epoched data?
CFG.binlister = PATHS.binlister(max(strfind(PATHS.binlister,filesep))+1:strfind(PATHS.binlister,'.txt')-1); % automatically get Binlister file name

%% perform checks before running script

disp( fullfile(PATHS.dataDir,PATHS.subjdatafile) )

% if directory folder doesn't exist, skip this subject
if ~exist(fullfile(PATHS.dataDir,PATHS.subjdatafile))
    return
end

disp('two')

if length(dir(fullfile(PATHS.dataDir,PATHS.subjdatafile))) ~= 1 % if more than one raw data found, error
    return
end

try
%% read data

disp('three')

clc; close all;
disp(['Processing ' PATHS.subjdatafile, '...']);
load(fullfile(PATHS.dataDir,PATHS.subjdatafile));

%% Create eventlist, assign bins, and extract epochs using ERPLAB

EEG.EVENTLIST = []; % clear EVENTLIST
% pop_squeezevents(EEG); % summarize events (ERPLAB)
% eeg_eventtypes(EEG) % summarize events (EEGLAB)

% create eventlist with ERPLAB and save it in results folder
outPath = fullfile('Events',[CFG.event0]); mkdir(outPath);
EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'},'Eventlist',fullfile(outPath,[EEG.subject '_eventlist_raw.txt']));

% assign bins with ERPLAB binlister
EEG = pop_binlister(EEG,'BDF',PATHS.binlister,'ExportEL',fullfile(outPath,[EEG.subject '_bins.txt']),'IndexEL',1,'SendEL2','EEG','UpdateEEG','on','Voutput','EEG');

% extract epochs with ERPLAB
if ~isempty(CFG.epochBaselineCorrect)
    EEG = pop_epochbin(EEG,CFG.epochRange,CFG.epochBaselineCorrect); % extract epochs with baseline correction
else
    EEG = pop_epochbin(EEG,CFG.epochRange,'none'); % no baseline correction
end

EEG.condition = {EEG.EVENTLIST.bdf.description}; % bin description
EEG = eeg_checkset(EEG);
EEG.comments = pop_comments(EEG.comments,'','Extracted epochs',1);
disp(EEG.comments(end,:))

% Check if time-locking events in epochs match: EEG.EVENTLIST.eventinfo.bini & EEG.EVENTLIST.eventinfo.bepoch
timelockeventsbepoch = [[EEG.EVENTLIST.eventinfo.bepoch] > 0];
if sum([EEG.EVENTLIST.eventinfo(timelockeventsbepoch).bini] == -1)
    error('Error! Check EEG.EVENTLIST.eventinfo.bini and EEG.EVENTLIST.eventinfo.bepoch!');
else
    EEG.comments = pop_comments(EEG.comments,'','Epoched data and events match!',1);
    disp(EEG.comments(end,:))
end

%% Artifact detection on epoched data with EEGLAB

disp('=========================')
disp(fullfile(PATHS.artifact,[CFG.event0 '_EEGreject'],[EEG.subject '.mat']))
disp('=========================')

if exist(fullfile(PATHS.artifact,[CFG.event0 '_EEGreject'],[EEG.subject '.mat']),'file') % reload artifact detection .mat file

    % load artifact detection and rejected epoch information from .mat file
    disp('Not rerunning artifact detection...');
    EEG.comments = pop_comments(EEG.comments,'','Overwriting existing EEG.reject with previously saved EEG.reject stored in .mat file',1);
    disp(EEG.comments(end,:))
    disp(fullfile([CFG.event0 '_EEGreject'],[EEG.subject '.mat']));
    load(fullfile(PATHS.artifact,[CFG.event0 '_EEGreject'],[EEG.subject '.mat']));
    EEG.reject = EEGRejectField;
    EEG = eeg_checkset(EEG);

else % run artifact detection

    if ~exist(fullfile(PATHS.artifact,[CFG.event0 '_EEGreject']), 'dir')
        mkdir(fullfile(PATHS.artifact,[CFG.event0 '_EEGreject']));
    end

    % Citation: Delorme, A., Sejnowski, T., & Makeig, S. (2007). Enhanced detection of artifacts in EEG data using higher-order statistics and independent component analysis. NeuroImage, 34(4), 1443-1449. doi:10.1016/j.neuroimage.2006.11.004
    disp('Detecting artifacts...');

    % exclude eye and EMG channels from artifact detection
    eyeEmgChans = {'SO1','IO1','SO2','IO2','LO1','LO2','CorsOut','CorsIns','CORRins','CORRout','ZYGup','ZYGlow','COORins','COORout'}; % chan to exclude
    allChannels = {EEG.chanlocs.labels}; % all channel labels
    channelsToDetectArtifact = find(~ismember(lower(allChannels),lower(eyeEmgChans))); % exclude channels above

    % reject by linear trend/variance (max slope, uV/epoch: 100; R-squred limit: 0.5)
    EEG = pop_rejtrend(EEG,1,channelsToDetectArtifact,EEG.pnts,100,0.5,0,0,0);
    disp(['Trials rejected via linear trend/variance: ' num2str(find(EEG.reject.rejconst))]);

    % reject by probability (5 SD)
    EEG = pop_jointprob(EEG,1,channelsToDetectArtifact,5,5,0,0,0,[],0);
    disp(['Trials rejected via probability (5 SD): ' num2str(find(EEG.reject.rejjp))]);

    % reject by spectra (detect muscle and eye movement)
    % muscle: -100 to 100 dB, 20 to 40 Hz
    % eye: -50 to 50 dB, 0 to 2 Hz
    EEG = pop_rejspec(EEG,1,'elecrange',channelsToDetectArtifact,'method','fft','threshold',[-50 50;-100 25],'freqlimits',[0 2;20 40],'eegplotcom','','eegplotplotallrej',0,'eegplotreject',0);
    disp(['Trials rejected via spectra: ' num2str(find(EEG.reject.rejfreq))]);

    % update EEG.reject.rejglobal and EEG.reject.rejglobalE fields with all rejected epochs
    EEG = eeg_rejsuperpose(EEG,1,1,1,1,1,1,1,1);

    % save EEG.reject field as a .mat file so don't need to rerun artifact detection again in the future
    EEGRejectField = EEG.reject;
    save(fullfile(PATHS.artifact,[CFG.event0 '_EEGreject'],[EEG.subject '.mat']),'EEGRejectField');
    EEG.comments = pop_comments(EEG.comments,'','Finished artifact detection and saved EEG.reject field',1);
    disp(EEG.comments(end,:))
end

% find([EEG.EVENTLIST.eventinfo.flag]);
EEG = pop_syncroartifacts(EEG,'Direction','eeglab2erplab'); % transfer EEG.reject (eeglab) info to EEG.EVENTLIST (erplab)
EEG = eeg_checkset(EEG);

%% Summarize artifact detection info

% save this subject's artifact information into artifacts directory
colNames = {'subject','epochs','rejLinear','rejProb','rejSpec','rejTotal','rejPerc','acceptTotal','acceptPerc','binlister'};
values = {EEG.subject,EEG.trials,length(find(EEG.reject.rejconst)),length(find(EEG.reject.rejjp)),length(find(EEG.reject.rejfreq)),length(find(EEG.reject.rejglobal)),round(length(find(EEG.reject.rejglobal))/EEG.trials*100,2),EEG.trials-length(find(EEG.reject.rejglobal)),round((EEG.trials-length(find(EEG.reject.rejglobal)))/EEG.trials*100,2),CFG.binlister};
T = cell2table(values,'VariableNames',colNames);

outPath = fullfile(PATHS.artifact,CFG.event0); mkdir(outPath);
writetable(T, fullfile(outPath,[EEG.subject '_bin_' CFG.binlister]),'delimiter',',');
EEG.comments = pop_comments(EEG.comments,'','Saved artifact information',1);
disp(EEG.comments(end,:))

% save/update all subjects' artifact info and update csv file
% gather all subject's artifact detection summary in a table as save as csv file in Binlister directory
% tempFiles = dir(fullfile(dataDirectoryAbsPath,'**','parameters','',['*_bin_artifactN_' binlister])); % find matching files recursively (only works for MATLAB 2016b onwards)
% gather all subject's artifact detection summary in a table as save as csv file in Binlister directory
artifactTableAllSubjects = table();
tempFiles = dir(fullfile(PATHS.artifact,CFG.event0,['*' CFG.binlister '*.txt'])); % find matching files recursively
filePaths = {tempFiles.name};
for i=1:length(filePaths)
    tempData = readtable(fullfile(outPath,filePaths{i}));
    artifactTableAllSubjects = [artifactTableAllSubjects; tempData];
end
artifactTableAllSubjects.reject = artifactTableAllSubjects.rejPerc > CFG.rejecttrialpercent; % if > trials rejected, mark subject to reject
artifactTableAllSubjects = [table([1:height(artifactTableAllSubjects)]') artifactTableAllSubjects]; % add subject number count
artifactTableAllSubjects.Properties.VariableNames{1} = 'subjectN';

writetable(artifactTableAllSubjects,fullfile(PATHS.artifact,CFG.event0,[filetoread_noExt,'_artifact_summary.csv']));
disp('Saved all subject''s artifact information to directory');

%% Save single-trial event information (to get design matrix)

el = EEG.EVENTLIST.eventinfo; % get ERPLAB eventlist structure
timeLockedEventsIdx = [el.bepoch] ~= 0; % find indices with time-locked events (bepoch is not 0)
allEvents = {el.binlabel}; % all event labels (B1, B2 etc.)
epochEvent = allEvents(timeLockedEventsIdx); % event/code for every epoch (including epochs with artifacts)
% find artifact and no-artifact epochs
artifactFlags = [el.flag];
cleanEpochs = find(~artifactFlags(timeLockedEventsIdx)); % epoch (indices) without artifacts

% see if my manual count of trials and EEG.trials match
if ((length(cleanEpochs) + sum(artifactFlags)) ~= EEG.trials) || (length(epochEvent) ~= EEG.trials) 
    error('Error! Trial numbers don''t match! Double check!');
end

% Save single-trial info (design matrix) as table
elT = struct2table(el);
elT.enable = []; elT.dura = []; elT.codelabel = []; % remove variables
% add subject id to table
C = cell(size(elT,1),1); % empty cell to store subject id
C(:) = {EEG.subject}; % fill empty cell with subject id
elT = [C elT]; % join cell with table
elT.Properties.VariableNames = {'subject','eventN','eventCode','binlabel','timeS','samplingPoint','artifactFlag','binIndicator','epochN'}; % rename variable names

% join design matrix intro if it exists
if exist(PATHS.designmat, 'dir')
    toread = dir(fullfile(PATHS.designmat,[EEG.subject '_designmat.csv']));
    if length(toread) == 1
        designmattable = readtable(fullfile(PATHS.designmat,toread.name),'TreatAsEmpty','NA');
        designmattable.subject = [];
        elT = outerjoin(elT,designmattable,'Type','left','Mergekeys',true);
    end
end

elT.bindescr = elT.binlabel; % assign bindescr to each epoch
for i=1:size(elT,1) %  for each event...
    bindescrTemp = elT{i,'bindescr'}{1}; % get binlabel
    if strcmpi(bindescrTemp,'""') % if empty, it's not time-locking event
        elT{i,'bindescr'}{1} = '';
    else  % if not empty, fill in bindescr from ERP
        elT{i,'bindescr'}{1} = EEG.condition{elT{i,'binIndicator'}};
    end
end

el_timelockevent = elT((elT.epochN ~= 0),:); % eventlist with just time-locked event info

% save all events as csv
outPath = fullfile('Events',[CFG.event0]); mkdir(outPath);
writetable(elT,fullfile(outPath,[EEG.subject '_event_all.csv'])) % all events

%% Identify epochs to exclude

% check if design matrix is okay...
if sum(el_timelockevent{:,'epochN'} < 0) || sum(el_timelockevent{:,'binIndicator'} < 0)
    error('Error! Check design matrix epoch number and bin number!');
end

% mark trials with artifact
el_timelockevent.include = el_timelockevent{:,'artifactFlag'} == 0;
% for each column/variable in design matrix, if any row has NaN, include = 0
for col=1:size(el_timelockevent,2)
    if isnumeric(el_timelockevent{:,col})
        el_timelockevent{find(isnan(el_timelockevent{:,col})),'include'} = 0;
    end
end

%% Save clean epochs and design matrix

% select only trials where include == 1
EEG = pop_select(EEG,'trial',find(el_timelockevent{:,'include'} == 1));
EEG.comments = pop_comments(EEG.comments,'','Select epochs where include == 1',1);
disp(EEG.comments(end,:))

designmatclean = el_timelockevent(el_timelockevent{:,'include'} == 1,:);
EEG.designmat = designmatclean;

% hausekeeping
EEG.PATHS = PATHS;
EEG.CFG_epoch = CFG;
EEG.history = [];
EEG.setname = [EEG.subject ' ' CFG.event0];

% save design matrix and EEG struct
outPath = fullfile(PATHS.epochdir,CFG.event0); mkdir(outPath);
writetable(designmatclean,fullfile(outPath,[EEG.subject '_designmat.csv'])) 

if CFG.epochSave
    if ~exist(outPath, 'dir')
        mkdir(outPath);
    end

    save(fullfile(outPath,[EEG.subject '.mat']), 'EEG');
end

%% Done

disp(['Finished subject ' filetoread_noExt, '!']);
dlmwrite(fullfile(PATHS.log,['Finished_Subject_' filetoread_noExt ' ' datestr(datetime('now'), tformat) '.txt']),['Finished!' ' ' datestr(datetime('now'))],'delimiter','');

catch ME % if errors when running script, catch them

disp(['Error with subject ' filetoread_noExt, '!']);
save(fullfile(PATHS.log,['Error_MException_' filetoread_noExt ' ' datestr(datetime('now'), tformat) '.mat']),'ME'); % save MException object

end

% clear
close all

end % end function
