function retrace(pulse, opts)
%RETRACE
% Re-trace the sub-sequence sampling with different margins left and right
% of center frame, as specified in the new FIT_OPT
%
% USAGE: pulse.retrace( fit_opts )

fits = [pulse.fits]; cells = [pulse.cells];

for i = 1:numel(fits)
    
    this_fit = fits(i);
    % Find corresponding cell
    this_cell = pulse.find_cells_with_fit(this_fit);
    
    num_frames = numel(nonans(this_cell.dev_time));
    newOpt = opts(this_fit.embryoID);
    center_frame = this_fit.center_frame;
    
    % Get margin frame
    [left_margin,pad_l] = max( ...
        [center_frame - newOpt.left_margin, 1] );
    [right_margin,pad_r] = min( ...
        [center_frame + newOpt.right_margin, num_frames] );
    this_fit.margin_frames = left_margin:right_margin;
    
    % Collect the fit curve
    x = this_cell.dev_time( left_margin:right_margin );
    fitted_y = lsq_gauss1d( ...
        [this_fit.amplitude, this_fit.center, this_fit.width], x);
    this_fit.raw = this_cell.(newOpt.to_fit)( left_margin: right_margin );
    this_fit.fit = fitted_y;
    this_fit.aligned_time = x - this_fit.center;
    
    % PAD
    if pad_l > 1
        fitted_y = [ensure_row(nan(1 - (center_frame - newOpt.left_margin), 1)), fitted_y];
        x = [ensure_row(nan(1 - (center_frame - newOpt.left_margin), 1)), x];
    end
    if pad_r > 1
        fitted_y = [fitted_y, nan(1, (center_frame + newOpt.right_margin) - num_frames + 1)];
        x = [x, nan(1, (center_frame + newOpt.right_margin) - num_frames)];
    end
    this_fit.aligned_time_padded = x - this_fit.center;
    this_fit.opt = newOpt;
    
end

end %retrace
