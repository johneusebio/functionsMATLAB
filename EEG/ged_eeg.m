function [results] = ged_eeg(cfg)
% ged_eeg
% Performs GED with S and R matrices
% 
% Type the function name without any inputs to get parameter information
% 
% Last modified by Hause Lin 12-05-19 5:13 PM hauselin@gmail.com

%% display function parameters if no arguments provided

if nargin == 0 % if no input provided, return function argument information
    disp('ged_eeg function parameters: ');
    struct('data','eeglab data struct',...
        'Rdata','eeglab data struct for R covariance matrix (default: uses S matrix eeglab struct)',...
        'singletrial','single-trial covariance (1 or 0; default: 1, yes)',...
        'Swin','S matrix time window',...
        'Rwin','R matrix time window',...
        'regularize','shrinkage/regularization (0 to 1; default 0)',...
        'components','components to compute and plot (default = 5)',...
        'rawdata','data X for computing component time series (c = wX)',...
        'plotfig','plot results figure number (default: 0)',...
        'plottime','plot xlim for time series (default: [min max])',...
        'verbose','print messages for debugging (default: 0)')
    return
end

%% set up default values

if ~isfield(cfg,'singletrial') 
    cfg.singletrial = 1; % default computes single-trial covariance
end

if ~isfield(cfg,'Rdata')
    cfg.Rdata = cfg.data; % default computes S and R covariance matrices using the same data
end

if ~isfield(cfg,'components')
    cfg.components = 5; % default precomputes and plots top 5 components
    cfg.timeseriescomponents = cfg.components;
else
    cfg.components = max(3,cfg.components); % fewest components to plot is 3
    cfg.timeseriescomponents = cfg.components;
end

if ~isfield(cfg,'plotfig')
    cfg.plotfig = 0; % default don't plot
end

if ~isfield(cfg,'regularize')
    cfg.regularize = 0; % default no regularization
end

if ~isfield(cfg,'verbose')
    cfg.verbose = 0; 
end

%% get time indices for R and S matrices

cfg.Rwinidx = dsearchn(cfg.Rdata.times',[cfg.Rwin(1) cfg.Rwin(end)]');
cfg.Swinidx = dsearchn(cfg.data.times',[cfg.Swin(1) cfg.Swin(end)]');

%% perform GED

cfg.data.data = double(cfg.data.data); % ensure double precision for eigendecomposition
cfg.Rdata.data = double(cfg.Rdata.data); % ensure double precision for eigendecomposition

if cfg.verbose
    disp(['Performing GED on ' num2str(cfg.data.trials) ' trials (S matrix)']);
    disp(['Performing GED on ' num2str(cfg.Rdata.trials) ' trials (R matrix)']);
end

% compute channel-by-channel covariance
[covR,covS] = deal(zeros(cfg.data.nbchan)); % initialize empty covariance matrices
if cfg.singletrial % 
    if cfg.verbose
        disp('Computing single-trial covariance matrices...');
    end
    
    if cfg.Rdata.trials == cfg.data.trials % if S and R matrices have same number of trials, only 1 loop needed
        for ti=1:cfg.data.trials % covariance matrix for each trial
            % R matrix
            tdat = squeeze(cfg.Rdata.data(:,cfg.Rwinidx(1):cfg.Rwinidx(2),ti)); 
            covR = covR + cov(tdat');

            % S matrix
            tdat = squeeze(cfg.data.data(:,cfg.Swinidx(1):cfg.Swinidx(2),ti)); 
            covS = covS + cov(tdat');
        end
    else % if S and R matrices don't' have same number of trials
        for ti=1:cfg.Rdata.trials % covariance matrix for each trial
            % R matrix
            tdat = squeeze(cfg.Rdata.data(:,cfg.Rwinidx(1):cfg.Rwinidx(2),ti)); 
            covR = covR + cov(tdat');
        end
        
        for ti=1:cfg.data.trials % covariance matrix for each trial
            % S matrix
            tdat = squeeze(cfg.data.data(:,cfg.Swinidx(1):cfg.Swinidx(2),ti)); 
            covS = covS + cov(tdat');
        end
    end
    
    covR = covR./cfg.Rdata.trials; % average covariance
    covS = covS./cfg.data.trials;
else % covariance matrix on avg data
    if cfg.verbose
        disp('Computing covariance matrices on averaged data (e.g., ERP)...');
    end
    % just 1 S and 1 R matrix!
    
    % R matrix
    avgdata = mean(cfg.Rdata.data,3); % compute averaged data (ERP)
    tdat = squeeze(avgdata(:,cfg.Rwinidx(1):cfg.Rwinidx(2))); 
    covR = cov(tdat');

    % S matrix
    avgdata = mean(cfg.data.data,3); % compute averaged data (ERP)
    tdat = squeeze(avgdata(:,cfg.Swinidx(1):cfg.Swinidx(2)));
    covS = cov(tdat');
end

% apply shrinkage regularization to R matrix
if cfg.regularize > 0
    covR = (1 - cfg.regularize) * covR + cfg.regularize * mean(eig(covR)) * eye(cfg.data.nbchan);
end

% perform eigendecomposition
[evecs,evals] = eig(covS,covR);
if ~isreal(evals) % if returns imaginary values, regularize shrink 1%
    disp('Regularizing by 1% because eigendecomposition returned imaginary values!');
    cfg.regularize = 0.01;
    covR = (1 - cfg.regularize) * covR + cfg.regularize * mean(eig(covR)) * eye(cfg.data.nbchan);
    [evecs,evals] = eig(covS,covR);
end

% sort eigenvectors
[evals,sidx] = sort(diag(evals),'descend');
evalsprop = evals ./ sum(evals) * 100; % proportion variance explained by each eigenvalue
evecs = evecs(:,sidx); % sort eigenvectors by eigenvalue

%% normalize evecs such that magnitude/vector norm is 1

% ensures component time series has similar amplitude as ERP 
evecs = evecs./sqrt(sum(evecs.^2,1));
% norm(evecs(:,1)) % should be 1

%% compute filter-forward model/activation pattern and flip sign

% compute activation pattern/topography for each eigenvector (a = wS, where S is covariance matrix S)
activationpatterns = [];
if cfg.components
    activationpatterns = zeros(cfg.data.nbchan,size(evecs,2));
    evecsignflip = zeros(1,size(evecs,2));
    for c=1:size(evecs,2)
        activationpatterns(:,c) = evecs(:,c)'*covS; % compute activation pattern/topography (a = wS, where S is covariance matrix S)
        % flip eigenvector sign if necessary
        [~,idx] = max(abs(activationpatterns(:,c)));  % find absolute max magnitude in activation pattern
        evecsignflip(c) = sign(activationpatterns(idx,c));
        evecs(:,c) = evecs(:,c) * evecsignflip(c); % possible sign flip
        activationpatterns(:,c) = activationpatterns(:,c)*evecsignflip(c); % possible sign flip
    end    
end
% note: we can also flip sign based on channel

%% compute component time series

% compute component time series (projections) (c = wX, where X is data matrix)
timeseriescomponents = [];
if cfg.timeseriescomponents
    if isfield(cfg,'rawdata')
        timeseriescomponents = evecs(:,1:cfg.timeseriescomponents)'*squeeze(mean(cfg.rawdata,3)); % comp_chan * chan_time
    else
        % apply weights to trial-averaged data
        timeseriescomponents = evecs(:,1:cfg.timeseriescomponents)'*squeeze(mean(cfg.data.data,3)); % comp_chan * chan_time
        % apply weights to single-trial data then average (same results as above)
        % dat2d = reshape(cfg.data.data,size(cfg.data.data,1),[]);
        % timeseriescomponents_singletrial = evecs(:,1:cfg.timeseriescomponents)' * dat2d;
        % timeseriescomponents_singletrial = reshape(timeseriescomponents_singletrial,cfg.timeseriescomponents,size(cfg.data.data,2),size(cfg.data.data,3));
        % timeseriescomponents_singletrial = mean(timeseriescomponents_singletrial,3);
    end    
end

%% visualize results

if cfg.plotfig
    
    figure(cfg.plotfig); clf
    
    % plot eigenspectrum
    subploteigspec = [1:cfg.components-2];
    subplot(3,cfg.components,subploteigspec)
    plot(1:length(evalsprop),evalsprop,'s-','linew',1,'markersize',5,'markerfacecolor','k')
    set(gca,'xlim',[1 length(evalsprop)])
    ylabel('Variance explained'); title(['Eigenvalues']);
    if cfg.singletrial
        t = '(single-trial covariance)';
    else
        t = '(ERP covariance)';
    end
    if length(evalsprop) > 30
        title({['Eigenvalues for top 30 (of ' num2str(length(evalsprop)) ' values)'] ['Shrink: ' num2str(cfg.regularize) ' ' t]},'fontsize',13);
        set(gca,'xlim',[1 30])
    else
        title({'Eigenvalues' ['Shrink: ' num2str(cfg.regularize) ' ' t]},'fontsize',13);
    end
    
    % plot input data topography
    subplot(3,cfg.components,subploteigspec(end)+1)
    tempdat = cfg.data.data(:,cfg.Swinidx(1):cfg.Swinidx(2),:);
    topoplotIndie(squeeze(mean(mean(tempdat,2),3)),cfg.data.chanlocs,'electrodes','labels');
    title({'GED S matrix data' [num2str(round(cfg.Swin(1))) '-' num2str(round(cfg.Swin(2)))]});
    colorbar; caxis([-max(abs(caxis)) max(abs(caxis))])
    
    if isfield(cfg,'rawdata')
        subplot(3,cfg.components,subploteigspec(end)+2)
        tempdat = cfg.rawdata(:,cfg.Swinidx(1):cfg.Swinidx(2),:);
        topoplotIndie(squeeze(mean(mean(tempdat,2),3)),cfg.data.chanlocs,'electrodes','on');
        title('Data for component time series');
        colorbar; caxis([-max(abs(caxis)) max(abs(caxis))])
    else
        subplot(3,cfg.components,subploteigspec(end)+2)
        tempdat = cfg.Rdata.data(:,cfg.Rwinidx(1):cfg.Rwinidx(2),:);
        topoplotIndie(squeeze(mean(mean(tempdat,2),3)),cfg.data.chanlocs,'electrodes','on');
        title({'GED R matrix data' [num2str(round(cfg.Rwin(1))) '-' num2str(round(cfg.Rwin(2)))]});
        colorbar; caxis([-max(abs(caxis)) max(abs(caxis))])
    end

    % plot activation patterns/topography (a = wS)
    for i=1:cfg.components
        subplot(3,cfg.components,cfg.components+i) % second row
        topoplotIndie(activationpatterns(:,i),cfg.data.chanlocs,'electrodes','on');
        title([ 'Component ' num2str(i) ])
        colorbar; caxis([-max(abs(caxis)) max(abs(caxis))])
        
    end
    
    % plot component time series (projections) (c = wX, where X is data matrix)
    for i=1:cfg.components
        subplot(3,cfg.components,2*cfg.components+i) % third row
        plot(cfg.data.times,timeseriescomponents(i,:));
        if isfield(cfg,'plottime')
            xlim([cfg.plottime(1) cfg.plottime(2)])
        else
            xlim([cfg.data.times(1) cfg.data.times(end)])
        end
        xlabel('Time')
        if i == 1 
            ylabel('Amplitude');
        end
    end
    try
        colormap viridis
    catch
        try
            colormap inferno
        catch
            colormap parula
        end
    end
end

%% return result

results = [];
results.R = covR;
results.S = covS;
results.evecs = evecs;
results.evals = evals;
results.evalsprop = evalsprop;
results.evecsignflip = evecsignflip;
results.activationpatterns = activationpatterns;
results.timeseriescomponents = timeseriescomponents;
results.times = cfg.data.times;
results.chanlocs = cfg.data.chanlocs;

cfgout = cfg;
cfgout = rmfield(cfgout,'data');
if isfield(cfg,'rawdata')
    cfgout = rmfield(cfgout,'rawdata');
end
if isfield(cfg,'Rdata')
    cfgout = rmfield(cfgout,'Rdata');
end
results.cfg = cfgout;

end

