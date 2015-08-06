function fits_output = find_fits_from_cell(pulse,cells)
%FIND_FITS_FROM_CELLS Return the Fits from the given CellObjs.
%
% USAGE: fits = pulses.find_fits_from_cells(cells);

embryoIDs = unique([cells.embryoID]);
assert( all(ismember(embryoIDs,[pulse.embryoID])),'Not all embryoIDs found in Pulse input');
assert( strcmpi(class(cells),'CellObj'),'Those aren''t cells...');

fits_output = cell(1,numel(embryoIDs));

for i = 1:numel(embryoIDs)

    cellsOI = cells([cells.embryoID] == embryoIDs(i));
    fits = [pulse.get_embryoID(embryoIDs(i)).fits];
    fits = fits( ismember([fits.cellID],[cellsOI.cellID]) );

    fits_output{i} = fits;
    
end

fits_output = [fits_output{:}];

end