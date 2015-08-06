function cells = get_cellID(pulse,cellID)
%GET_CELLID Return the cellObj(s) with the given cellID(s) from a single pulse.
% 
% USAGE: cells = pulse.get_cellID(cellID);

assert(numel(pulse) == 1)

cells = pulse.cells;
cells = cells( ismember([cells.cellID],cellID) );

end