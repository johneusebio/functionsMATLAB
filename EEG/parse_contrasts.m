function [contrasts, contrasts_labels] = parse_contrasts(filename, labels)
%parse_contrasts - Description
%
% Syntax: contrasts = parse_contrasts(input_txt)
%
% Long description

%% Perform checks
if isa(filename, 'char')
    if ~exist(filename, 'file')
        disp(['"' filename '" is not a valid file path']);
        return
    end
else
    disp('The "filename" must be in "char" format.');
    return
end

if nargin > 1
    if ~isa(labels, 'cell')
        disp('"labels" must be a cell of char variables');
        return
    end
end

input_txt = importdata(filename);

contrasts = struct();

for line = 1:size(input_txt,1)
    txt_line = input_txt{line,1};

    % check if that particular contrast has been commented out
    if startsWith(txt_line, '%')
        continue
    end

    txt_line = strsplit(txt_line, ':');

    contrast_nm   = txt_line{1};
    contrast_bins = strsplit(strtrim(txt_line{2}), ' ');

    neg_Ind = find(strcmp(contrast_bins, '-'));
    bin_Ind = ~contains(contrast_bins, {'+', '-'});

    contrast_bins(neg_Ind+1) = strcat(contrast_bins(neg_Ind), contrast_bins(neg_Ind+1));
    contrast_bins = contrast_bins(bin_Ind);
    contrast_bins = cellfun(@str2num, contrast_bins);

    contrast_weights = -1 * (contrast_bins < 0);
    contrast_bins    = abs(contrast_bins);

    contrast_weights(contrast_weights==0) = 1;

    contrasts(line).name    = contrast_nm;
    contrasts(line).bins    = contrast_bins;
    contrasts(line).weights = contrast_weights;
end

if nargin > 1
    contrasts_labels        = struct();
    contrasts_labels.bins   = unique([contrasts(:).bins]);
    contrasts_labels.labels = unique(labels([contrasts(:).bins]));
end

end