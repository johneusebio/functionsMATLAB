function bin_count_table = bin_counter(bins, EEG, count_thr)
%bin_counter - Description
%
% Syntax: nbins = bin_counter(bins, ERP)
%
% Long description

bin_count       = sum(bins == EEG.designmat.binIndicator, 1);
bin_count_table = table(bins', bin_count', 'VariableNames', {'bins', 'count'});

if nargin > 2
    bin_count_low = bin_count < count_thr;

    if ~exist('Warnings', 'dir')
        mkdir('Warnings');
    end

    if any(bin_count_low)
        fileID = fullfile('Warnings', ['WARNING_nBins_' EEG.subject '.txt']);
        writetable(bin_count_table, fileID);
    end
end

end