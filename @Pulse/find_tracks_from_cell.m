function tracks_output = find_tracks_from_cell(pulse,cells)
%FIND_FITS_FROM_CELLS Return the Tracks from the given CellObjs.
%
% USAGE: fits = pulses.find_tracks_from_cell(cells)

embryoIDs = unique([cells.embryoID]);
assert( all(ismember(embryoIDs,[pulse.embryoID])),'Not all embryoIDs found in Pulse input');
assert( strcmpi(class(cells),'CellObj'),'Those aren''t cells...');

tracks_output = cell(1,numel(embryoIDs));

for i = 1:numel(embryoIDs)

    cellsOI = cells([cells.embryoID] == embryoIDs(i));
    tracks = [pulse.get_embryoID(embryoIDs(i)).tracks];
    tracks = tracks( ismember([tracks.cellID],[cellsOI.cellID]) );

    tracks_output{i} = tracks;
    
end

tracks_output = [tracks_output{:}];

end