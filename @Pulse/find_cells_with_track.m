function cells_output = find_cells_with_track(pulse,tracks)
%FIND_CELL_WITH_FIT Return the CellObjs that contains the Track objects.
%
% USAGE: cells = pulses.find_cell_with_track(tracks);

embryoIDs = unique([tracks.embryoID]);
assert( all(ismember(embryoIDs,[pulse.embryoID])),'Not all embryoIDs found in Pulse input');
assert( strcmpi(class(tracks),'Track'),'Those aren''t tracks...');

cells_output = cell(1,numel(embryoIDs));

for i = 1:numel(embryoIDs)

    fitsOI = tracks([tracks.embryoID] == embryoIDs(i));
    cells = [pulse.get_embryoID(embryoIDs(i)).cells];
    cells = cells( ismember([cells.cellID],[fitsOI.cellID]) );

    cells_output{i} = cells;
    
end

cells_output = [cells_output{:}];

end