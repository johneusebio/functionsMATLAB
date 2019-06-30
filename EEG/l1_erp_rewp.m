function [ERP] = l1_erp_rewp(subject_file, contrast_file)

clc

%% get the subject file name without the extension
[~, filetoread_noExt, ~] = fileparts(subject_file);

% subj = '35';
PATHS.cwd        = pwd;
PATHS.datadir    = fullfile(PATHS.cwd,'Epochs','stimulus');
PATHS.datafile   = fullfile(PATHS.datadir,subject_file);
PATHS.resultsdir = fullfile(PATHS.cwd,'Results');

if ~exist(PATHS.resultsdir, 'dir')
    mkdir(PATHS.resultsdir);
end

%% perform checks

% if ~exist(PATHS.datafile, 'file')
%     return
% end

%% load data and select channels

disp(['Processing subject ' filetoread_noExt])

load(PATHS.datafile);
EEG = pop_select(EEG,'nochannel',{'M1','M2','SO2','IO2'});
EEG.data = double(EEG.data);

% figure(1); clf
% topoplotIndie(EEG.data(:,400,200),EEG.chanlocs,'electrodes','labels');

%% low pass 

lowpass = 30;
if lowpass
    EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Boundary','boundary','Cutoff',lowpass,'Design','butter','Filter','lowpass','Order',2,'RemoveDC','off');
end

%% load contrasts
[contrasts, contrasts_labels] = parse_contrasts(contrast_file, EEG.condition);

%% compute ERPs

ERP = struct();

% create the ERP data for the uncontrasted bins
for bb = 1:length((contrasts_labels.labels))
    bb_label = contrasts_labels.labels{bb};

    if bb == 1
        cfg = []; cfg.data = EEG; cfg.designmat = EEG.designmat;
        cfg.var = 'bindescr'; cfg.filter = {bb_label};
        idx = get_epochidx(cfg);
        ERP.data = mean(EEG.data(:,:,idx.idx),3);
        ERP.condition = {cfg.filter};
        ERP.erpinfo = {idx};
    else
        cfg = []; cfg.data = EEG; cfg.designmat = EEG.designmat;
        cfg.var = 'bindescr'; cfg.filter = {bb_label};
        idx = get_epochidx(cfg);
        ERP.data(:,:,end+1) = mean(EEG.data(:,:,idx.idx),3);
        ERP.condition{end+1} = cfg.filter;
        ERP.erpinfo{end+1} = idx;
    end
end

% add on ERP data for the specified contrasts
for cc = 1:length(contrasts)

    % check if the participant has enough bins to warrant inclusion
    bin_count_table = bin_counter(contrasts(cc).bins, EEG, 10);

    if any(bin_count_table.count == 0)
        disp(['WARNING: Skipping contrast "' contrasts(cc).name '" due to insufficient number of bins.']);
        disp('Please check the corresponding file in the Warnings directory for details.');
        continue
    end

    bin_Ind = [];
    for bb = contrasts(cc).bins
        bin_Ind = [bin_Ind find(ismember(contrasts_labels.bins, bb))];
    end

    disp(bin_Ind)
    % find(ismember([contrasts(cc).bins], [contrasts_labels.bins]));

    contrast_matrix = ones(size(ERP.data(:,:, bin_Ind)));

    for ww = 1:length(contrasts(cc).weights)
        contrast_matrix(:,:,ww) = contrasts(cc).weights(ww);
    end

    ERP.data(:,:, end+1) = sum(ERP.data(:,:,bin_Ind) .* contrast_matrix, 3);
    ERP.condition{end+1} = {contrasts(cc).name};
    ERP.erpinfo{  end+1} = idx;
end

% include simple average
ERP.data(:,:,end+1) = mean(EEG.data,3);
ERP.condition{end+1} = {'avgsimple'};
ERP.erpinfo{end+1} = idx;

% copy the EEG fields to the ERP structure
ERP = copyEEGfields(EEG,ERP);

% save
savefile(fullfile(PATHS.resultsdir),[filetoread_noExt '.mat'],ERP);

%% EEG minus ERP

% figure(1); clf
% plot(ERP.times,ERP.data(get_chanidx(ERP,'Pz'),:,3),'linew',3)
% hold on
% EEGminusERP = EEG.data - squeeze(ERP.data(:,:,3));
% % EEGminusERP_avg = mean(EEGminusERP,3);

% trial = randi(200);
% plot(ERP.times,squeeze(EEG.data(get_chanidx(ERP,'Pz'),:,trial)),'linew',3);
% plot(ERP.times,squeeze(EEGminusERP(get_chanidx(ERP,'Pz'),:,trial)),'linew',3);
% xlim([-200 800]); legend({'erp' 'trial' 'trial minus erp'})

%% plot ERPs and topography

% help pop_plottopo
% pop_plottopo(ERP,[1:ERP.nbchan],'abc',1)
% topoplotIndie(ERP.data(:,dsearchn(ERP.times',[200]'),1),ERP.chanlocs,'electrodes','labels');

%% clean up

disp(['Finish subject ' filetoread_noExt])
close all

end % function end