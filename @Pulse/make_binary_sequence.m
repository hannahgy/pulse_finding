function binary = make_binary_sequence(pulse)
%MAKE_BINARY_SEQUENCE Uses width_frames to generate a binary
% sequence of pulses
% USAGE: binary_seq = pulse.make_binary_sequence;

if numel(pulse) > 1
    error('Only 1 Pulse at a time.');
end

cells = [pulse.cells];
fits = [pulse.fits];

% Preallocate
binary = zeros( numel(cells(1).dev_time), max( [cells.cellID] ));
% Filter relevant fits
fits = fits.get_fitID( [cells.fitID] );

for i = 1:numel(cells)
    
    this_cell_fits = fits.get_fitID( cells(i).fitID );
    if ~isempty(this_cell_fits)
        for j = 1:numel( this_cell_fits )
            binary( this_cell_fits(j).width_frames, cells(i).cellID ) = ...
                binary( this_cell_fits(j).width_frames, cells(i).cellID ) + 1;
            %                             this_cell_fits(j).cluster_label;
        end
    end
    
end

end % make_binary_sequence