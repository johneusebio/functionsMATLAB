function [] = trim_EGG(EEG_path, output, start_bin, end_bin)
%trim_restEGG - Description
%
% Syntax: [] = trim_restEGG(EEG, start_bin, stop_bin)
%
% Long description

if exist(EEG_path, 'file')
    disp('EEG must be a valid filepath.');
    return
end

load(EEG_path);

%% identify the latencies where the specified bins occur
event_latencies = [EEG.event.latency]
event_types     = cellfun(@str2num, {EEG.event.type}, 'UniformOutput', false);
event_types     = [event_types{:}];

start_latency   = event_latencies(event_types == start_bin) - 1;
end_latency     = event_latencies(event_types == end_bin) + 1;

keep_range = [1:start_latency, end_latency:size(EEG.times,2)];

%% begin removing times outsite of the specified range
EEG.times(keep_range)   = [];
EEG.data (:,keep_range) = [];

%% save to the output file path
save(output, 'EEG');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hdr    = ft_read_header(EEG ); 
% data   = ft_read_data(EEG, 'header', hdr );
% events = ft_read_event(EEG, 'header', hdr );