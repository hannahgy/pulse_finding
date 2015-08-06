function first_fits = get_first_fit(pulse)
%GET_FIRST_FIT Retrieve the first pulse (in time) from a cell
%
% SYNOPSIS: first_fits = cells.get_first_fit(fits);

cells = [pulse.cells];
% fits = [pulse.fits];

num_cells = numel(cells);
first_fits = cell(1,num_cells);
for i = 1:num_cells
    
    this_cell = cells(i);
    this_fits = pulse.get_embryoID( ...
        this_cell.embryoID ).find_fits_from_cell( this_cell );
%     this_fits = fits.get_stackID( this_cell.stackID );
%     this_fits = fits([fits.stackID] == this_cell.stackID);
    this_fits = this_fits.sort('center');
    if numel(this_fits) > 0
        first_fits{i} = this_fits(1);
    end
    
end

end % get_first_fit