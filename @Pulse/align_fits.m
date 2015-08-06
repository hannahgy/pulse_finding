function align_fits(pulse,name)
%ALIGN_PEAKS Aligns the global maxima of a given array of
%FITTED objects
% Will return also a given measurement aligned according to the
% maxima. Updates the FITTED structures contained in Pulse.
%
% SYNOPSIS: pulse.align_fits(name);
%           pulse.align_fits(name,measurement);
%

fits = [pulse.fits];
cells = [pulse.cells];

num_fits = numel(fits);
durations = cellfun(@numel, {fits.margin_frames} );

switch name
    case 'area'
        cell_name = 'area_sm';
    case 'myosin'
        cell_name = 'myosin_intensity';
    otherwise
        cell_name = name;
end

% Check whether measurement is numeric or cell array
cellArrFlag = iscell(cells(1).(cell_name));
vectorArrFlag = isscalar(cells(1).(cell_name));

for i = 1:num_fits
    
    this_fit = fits(i);
    
    if vectorArrFlag
        m = measurement(i);
    else
        
        frames = this_fit.margin_frames;
        opt = this_fit.opt;
        l = opt.left_margin; r = opt.right_margin;
        % center_idx indicates the index of the aligned matrix, not
        % the frame corresponding to the Gaussian center
        center_idx = l + 1;
        
        % If there is a tie, returns the first
        [~,max_idx] = max(this_fit.fit);
        left_len = max_idx - 1;
        
        if ~cellArrFlag
            m = nan(1, l + r + 1); % Make the padded vector
        else
            m = cell(1, l + r + 1);
        end
        
        lb = center_idx - left_len;
        ub = min(center_idx - left_len + durations(i) - 1, max(durations) );
        
        if numel(frames) - numel(lb:ub) == 1,
            frames = frames(1:end-1);
        end
        
        
        m( lb: ub) = ensure_row( ...
            pulse.find_cells_with_fit(this_fit).(cell_name)( frames ));
        
    end
    
    fits(i).(name) = m;
    
end

end % align_fits