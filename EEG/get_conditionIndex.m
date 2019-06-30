function condIndex = get_conditionIndex(filepath, conditions, delimiter)
%condIndex - Description
%
% Syntax: condIndex = get_conditionIndex(filepath, delimiter, conditions)
% 
% This function identifies which condition a given filepath belongs to. 
% The filepath must be entered as a char. It can be either a single file or a 
%   full filepath. It does not matter, so long as it's a char.
% The conditions must be entered either as a char or an array of multiple char.
% 
% The output will be an integer indicating which entry in the conditions array 
% the filename matches. If there is no match, 0 is returned.

if ~isa(filepath, 'char')
    disp('The specified filepath must be a char.');
    return
end

if ~isa(delimiter, 'char')
    disp('The specified delimiter must be a char.');
    return
end

if isa(conditions, 'char')
    conditions = {conditions};
elseif ~isa(conditions, 'cell')
    disp('Specified conditions must be entered as either a char or a cell.');
    % return
elseif isa(conditions, 'cell')
    condition_classes =  unique(cellfun(@class, conditions, 'UniformOutput', false));

    if length(condition_classes) > 1
        disp('The specified condition contain data of more than one class. All must be char.')
        % return
    elseif ~isa(condition_classes{1}, 'char')
        disp('The specified conditions are not in char format.')
        % return
    end
end

[~, filename, ~] = fileparts(filepath);

filename_chunks  = split(filename, delimiter);
condition_member = ismember(conditions, filename_chunks);

if sum(condition_member) > 0
    condIndex = find(condition_member);

    if length(condIndex) > 1
        disp(['More than one condition was found in the file name: ', filename])
        return
    end
else
    condIndex = 0;
end

end