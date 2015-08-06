function f = find_non_edge(pulse)
%FIND_NON_EDGE Finds the fits that are not on the edge of the
% segmented embryo.
%
% USAGE: fits_non_edge = pulses.find_non_edge;

f = [];
for i = 1:numel(pulse);
    
    fitsOI = pulse(i).fits;
    cellsOI = pulse(i).cells;
    cIDs = [cellsOI.cellID];
    
    tref = pulse(i).input.tref;
    
    neighborhood = cat(2,cellsOI.identity_of_neighbors_all);
    neighborhood = neighborhood(tref,:);
    
    N = cellfun(@(x) numel(x(ismember(x,cIDs))), neighborhood);
    
    foo = pulse(i).find_fits_from_cell(cellsOI(N > 4));
    f = [f foo];
    
end

end % find_non_edge
