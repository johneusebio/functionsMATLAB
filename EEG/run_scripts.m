clear, close all, %clc;
addpath('C:\Users\enter\OneDrive - University of Toronto\software\eeg_lab')
addpath('C:\Users\enter\OneDrive - University of Toronto\software\fieldtrip-20190611')

% M = 3;                     % M specifies maximum number of workers for parallel processing

%% load and close eeglab
eeglab
figHandles = get(groot, 'Children');
figHandles = {figHandles.Name};
for ii = strfind(figHandles, 'EEGLAB')
    close(figHandles{ii{1}});
end

data.path    = 'C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\data\lvl1\RawData';
data.ids     = strsplit(num2str(1001:1059,'%01d '),' '); % generate character subject ids
data.exclude = {'1035'}; % subjects which are to be excluded from the analysis
data.use     = data.ids(~ismember(data.ids, data.exclude)); % use only subjects who aren't excluded

%% define condition data for the experiment
data.conds           = struct();
data.conds.name      = {'rsvp', 'tap', 'audio'};
data.conds.binlist   = {'C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\experiment\analyze_data\scripts\support_files\rsvp_binlist.txt', ...
                        'C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\experiment\analyze_data\scripts\support_files\tap_binlist.txt' , ...
                        'C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\experiment\analyze_data\scripts\support_files\audio_binlist.txt'};
data.conds.contrasts = {'C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\experiment\analyze_data\scripts\support_files\rsvp_contrast.txt', ...
                        [], []};
data.conds.preproc_cfg = {[], [], ...
                          'C:\Users\enter\OneDrive - University of Toronto\Labs\Inzlicht Lab\outside_project\experiment\analyze_data\scripts\support_files\audio_preproc_cfg.txt'};

for subj_n=8:length(data.use)
    subj = data.use(subj_n);

    %% cnt_to_mat (convert cnt files to mat files)
    cfg = struct();
    cfg.rawdatadir = fullfile(data.path, subj{1}, 'eeg');
    cfg.matdatadir = fullfile(data.path, subj{1}, 'eegToMat');

    if ~exist(cfg.matdatadir, 'dir')
        mat_files = dir(cfg.matdatadir);
        mat_files = {mat_files.name};
        mat_files = mat_files(~cellfun(@isempty, strfind(mat_files, '.mat')));
        [~, mat_files, ~] = cellfun(@fileparts, mat_files, 'UniformOutput', false);

        raw_files = dir(cfg.rawdatadir);
        raw_files = {raw_files.name};
        raw_files = raw_files(~cellfun(@isempty, strfind(raw_files, '.cnt')));
        [~, raw_files, ~] = cellfun(@fileparts, raw_files, 'UniformOutput', false);

        cnt_to_mat(cfg);
    end

    eegToMat_list = dir(fullfile(cfg.matdatadir, '*.mat'));
    eegToMat_list = {eegToMat_list.name};

    for ss_file=eegToMat_list

        %% find which condition the current file belongs to, and grab the appropriate binlist
        ss_condInd = get_conditionIndex(ss_file{1}, data.conds.name, '_');
        
        if ss_condInd == 0
            continue
        end
        
        % s1_preprocessEEG
        s1_preprocessEEG(ss_file{1}, cfg, data.conds.preproc_cfg{ss_condInd});
        
        % then use s2_checkICA.m script to identify and remove IC components for each subject
        s2_checkICA(cfg, ss_file{1});
        
        %% l0_epoch_feedback.m
        ss_binlist   = data.conds.binlist{ss_condInd};
        l0_epoch_feedback(ss_file{1}, cfg, ss_binlist);
        
        %% l1_erp_rewp
        ss_contrasts = data.conds.contrasts{ss_condInd};
        if ~isempty(ss_contrasts)
            l1_erp_rewp(ss_file{1}, ss_contrasts);
        end
    end
end

%% John: Hause said the following scripts are for multivariate analyses and I won't need them at this stage.

%% l1ged_singletrial_vs_avg_covariance

% exclude = {'29'};
% subjectids = strsplit(num2str(1:46,'%01d '),' '); % generate character subject ids
% for s=subjectids
%     if ~ismember(exclude,s)
%         l1ged_singletrial_vs_avg_covariance(s{1});
%     end
% end

% %% l1ged_rewp_singletrial_vs_avg_covariance

% exclude = {'29'};
% subjectids = strsplit(num2str(1:46,'%01d '),' '); % generate character subject ids
% for s=subjectids
%     if ~ismember(exclude,s)
%         l1ged_rewp_singletrial_vs_avg_covariance(s{1});
%     end
% end

% %% end