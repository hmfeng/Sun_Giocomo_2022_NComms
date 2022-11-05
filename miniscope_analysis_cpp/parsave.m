function parsave(fpath, fname, PlaceCells_both, Tinfo)
% this function is used for save files within a parfor loop
if ~exist('Tinfo', 'var') || isempty(Tinfo) 
    save(fullfile(fpath, fname), 'PlaceCells_both');
else
    save(fullfile(fpath, fname), 'PlaceCells_both', 'Tinfo');
end
end