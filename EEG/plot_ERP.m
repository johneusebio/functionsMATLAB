function plot_ERP(subjFile, conds)
% plotERP - Description
%
% Syntax:plotERP(subjFile, conds)
%
% Long description

if ~isa(subjFile, 'char')
    disp('"subjFile" must be a char.');
    return
end

disp(['Loading ', subjFile, '...']);
load(fullfile(pwd, 'Results', subjFile));

if nargin > 1
    if ~isa(conds, 'double')
        disp('"conds" much be a double.');
        return
    end

    ERP.data = ERP.data(:,:,conds);
    ERP.condition = ERP.condition(conds);
    ERP.erpinfo = ERP.erpinfo(conds);
end

pop_plottopo(ERP,[1:ERP.nbchan], insertBefore(subjFile, '_', '\'),1);
topoplotIndie(ERP.data(:,dsearchn(ERP.times',[200]'),1),ERP.chanlocs);

end