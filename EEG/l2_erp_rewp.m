clear, clc, close all
cd '/Users/hause/Dropbox/Working Projects/Akina Physical Effort/Zaibas RewP study/Scripts/'

%% load subject erp data

PATHS = struct();
PATHS.cwd = pwd;
PATHS.resultsdir = fullfile(PATHS.cwd,'Results');
PATHS.datadir = fullfile(PATHS.resultsdir,'subj_erp_reward_noreward');
PATHS.data = dir(fullfile(PATHS.datadir,'*.mat'));

% load subject data
s=1;
tempdat = load(fullfile(PATHS.datadir,PATHS.data(s).name));
ERP = tempdat.ERP;
ERP.subj = {ERP.subject};
ERP.datasubj = ERP.data;
for s=2:length(PATHS.data)
    tempdat = load(fullfile(PATHS.datadir,PATHS.data(s).name));
    tempdat = tempdat.ERP;
    ERP.datasubj(:,:,:,s) = tempdat.data;
    ERP.subj(s) = {tempdat.subject};
end

% compute grand avg erps
ERP.data = nanmean(ERP.datasubj,4);
% ERP.datasub = [];

%% plot rewp erp and topography

chan2plot = 'FCz';
chan2plotidx = get_chanidx(ERP,'FCz');
figure(1); clf
subplot(2,2,1)
plot(ERP.times,squeeze(ERP.data(chan2plotidx,:,1)),...
    ERP.times,squeeze(ERP.data(chan2plotidx,:,2)),...
    ERP.times,squeeze(ERP.data(chan2plotidx,:,3)),...
    ERP.times,squeeze(ERP.data(chan2plotidx,:,4)),...
    'linew',2');
xlim([-200 800]);
xlabel('Time (s)','fontsize',16); ylabel('Amplitude','fontsize',16); title(chan2plot,'fontsize',18);
legend({'reward' 'noreward' 'reward-noreward' 'noreward-reward'},'fontsize',12,'location','northwest');

erpComponentWins = struct();
erpComponentWins.rewp = get_peakinfo(squeeze(ERP.data(chan2plotidx,:,3)),ERP.times,[150 350],'max',50,true);
erpComponentWins.rewp.chan = chan2plot;

figure(1)
y = get(gca,'ylim');
rectangle('Position',[erpComponentWins.rewp.peakstart y(1) erpComponentWins.rewp.duration diff(y)],'FaceColor',[0.5 0.5 0.5 0.3],'linestyle','none');

figure(1)
% plot topography
for t=1:3
    subplot(2,2,t+1)
    temptimes = [erpComponentWins.rewp.peakstart erpComponentWins.rewp.peakend];
    timeidx = dsearchn(ERP.times',temptimes');
    tempdat = squeeze(mean(ERP.data(:,timeidx(1):timeidx(2),t),2)); 
    topoplotIndie(tempdat,ERP.chanlocs);
    colorbar; cbsymmetric();
    title({ERP.condition{t}{1} [num2str(round(temptimes(1))) '-' num2str(round(temptimes(2))) ' ms']},'fontsize',16);
end

maximizefig()
save_figure(PATHS.resultsdir,'erpsTopo_rewp.png');

%% plot single-subject ERPS and topography

figure(2); clf
for s=1:size(ERP.datasubj,4)
    subplotgrid(size(ERP.datasubj,4),s);
    plot(ERP.times,squeeze(ERP.datasubj(chan2plotidx,:,1,s)),...
    ERP.times,squeeze(ERP.datasubj(chan2plotidx,:,2,s)),...
    ERP.times,squeeze(ERP.datasubj(chan2plotidx,:,3,s)),...
    'linew',2');
    xlim([-200 800]);
    y = get(gca,'ylim');
    rectangle('Position',[erpComponentWins.rewp.peakstart y(1) erpComponentWins.rewp.duration diff(y)],'FaceColor',[0.5 0.5 0.5 0.3],'linestyle','none');
    title([ERP.subj{s} ' ' chan2plot],'fontsize',11)
    if s == 1
        xlabel('Time (s)','fontsize',10); ylabel('Amplitude','fontsize',10);
    end    
end

maximizefig()
save_figure(PATHS.resultsdir,'erps_rewp_subj.png');

figure(3); clf
temptimes = [erpComponentWins.rewp.peakstart erpComponentWins.rewp.peakend];
timeidx = dsearchn(ERP.times',temptimes');
for s=1:size(ERP.datasubj,4)
    subplotgrid(size(ERP.datasubj,4),s);    
    tempdat = squeeze(mean(ERP.datasubj(:,timeidx(1):timeidx(2),3,s),2)); 
    topoplotIndie(tempdat,ERP.chanlocs);
    colorbar; cbsymmetric();
    title({[ERP.subj{s} ' ' ERP.condition{3}{1}] [num2str(round(temptimes(1))) '-' num2str(round(temptimes(2))) ' ms']},'fontsize',10);
end

maximizefig()
save_figure(PATHS.resultsdir,'topo_rewp_subj.png');

%% plot single-subject ERPS and topography

figure(2); clf
for s=1:size(ERP.datasubj,4)
    subplotgrid(size(ERP.datasubj,4),s);
    plot(ERP.times,squeeze(ERP.datasubj(chan2plotidx,:,1,s)),...
    ERP.times,squeeze(ERP.datasubj(chan2plotidx,:,2,s)),...
    ERP.times,squeeze(ERP.datasubj(chan2plotidx,:,3,s)),...
    'linew',2');
    xlim([-200 800]);
    y = get(gca,'ylim');
    rectangle('Position',[erpComponentWins.rewp.peakstart y(1) erpComponentWins.rewp.duration diff(y)],'FaceColor',[0.5 0.5 0.5 0.3],'linestyle','none');
    title([ERP.subj{s} ' ' chan2plot],'fontsize',11)
    if s == 1
        xlabel('Time (s)','fontsize',10); ylabel('Amplitude','fontsize',10);
    end    
end

maximizefig()
save_figure(PATHS.resultsdir,'erps_rewp_subj.png');

figure(3); clf
temptimes = [erpComponentWins.rewp.peakstart erpComponentWins.rewp.peakend];
timeidx = dsearchn(ERP.times',temptimes');
for s=1:size(ERP.datasubj,4)
    subplotgrid(size(ERP.datasubj,4),s);    
    tempdat = squeeze(mean(ERP.datasubj(:,timeidx(1):timeidx(2),3,s),2)); 
    topoplotIndie(tempdat,ERP.chanlocs);
    colorbar; cbsymmetric();
    title({[ERP.subj{s} ' ' ERP.condition{3}{1}] [num2str(round(temptimes(1))) '-' num2str(round(temptimes(2))) ' ms']},'fontsize',10);
end

maximizefig()
save_figure(PATHS.resultsdir,'topo_rewp_subj.png');

%% plot grdavg erp and topography

chan2plot = {'FCz' 'Cz' 'Pz'};
figure(4); clf
subplot(2,2,1)
for c=chan2plot
    chan2plotidx = get_chanidx(ERP,c{1});
    hold on
    plot(ERP.times,squeeze(mean(ERP.data(chan2plotidx,:,1:2),3)),'linew',2');
end
xlim([-200 800]);
xlabel('Time (s)','fontsize',16); ylabel('Amplitude','fontsize',16); title('grand avg','fontsize',18);
legend(chan2plot,'fontsize',12,'location','northwest');

tempdat = mean(mean(ERP.data(get_chanidx(ERP,chan2plot),:,1:2),3),1);
erpComponentWins.grdavg = get_peakinfo(tempdat,ERP.times,[100 500],'max',200,false);
erpComponentWins.grdavg.chan = chan2plot;

% manual window
erpComponentWins.grdavg.peakstart = 275;
erpComponentWins.grdavg.peakend = 425;
erpComponentWins.grdavg.duration = erpComponentWins.grdavg.peakend - erpComponentWins.grdavg.peakstart;

figure(4)
y = get(gca,'ylim');
rectangle('Position',[erpComponentWins.grdavg.peakstart y(1) erpComponentWins.grdavg.duration diff(y)],'FaceColor',[0.5 0.5 0.5 0.3],'linestyle','none');

figure(4)
% plot topography
for t=1:3
    subplot(2,2,t+1)
    temptimes = [erpComponentWins.grdavg.peakstart erpComponentWins.grdavg.peakend];
    timeidx = dsearchn(ERP.times',temptimes');
    tempdat = squeeze(mean(ERP.data(:,timeidx(1):timeidx(2),t),2)); 
    topoplotIndie(tempdat,ERP.chanlocs);
    colorbar; cbsymmetric();
    title({ERP.condition{t}{1} [num2str(round(temptimes(1))) '-' num2str(round(temptimes(2))) ' ms']},'fontsize',16);
end

maximizefig()
save_figure(PATHS.resultsdir,'erpsTopo_grdavg.png');

%% save erp component info

savefile([PATHS.resultsdir],'l2_erp_component_wins.mat',erpComponentWins);

%% end
