function fit = fit_gaussians(cells,opts)
%FIT_GAUSSIANS Fit multiple-Gaussians with a F-test stop
% Uses LSQCURVEFIT

num_embryos = numel(unique([cells.embryoID]));
num_cells = hist([cells.embryoID],1:num_embryos);

if numel(opts) ~= num_embryos
    error('Each embryo must have an associated fitting option structure.');
end

[cells(1:sum(num_cells)).flag_fitted] = deal(0);
num_frames = numel(cells(1).dev_time);

fitID = 0; last_embryoID = cells(1).embryoID; fit = [];

for i = 1:sum(num_cells)
    
    this_cell = cells(i);
    stackID = this_cell.stackID;
    if this_cell.embryoID ~= last_embryoID, fitID = 0; end % reset fitID count
    
    % Get relevant indices
    embryoID = this_cell.embryoID;
    % Get specific options
    opt = opts(1);
    
    Y = this_cell.(opt.to_fit); % curve to be fit
    t = this_cell.dev_time;     % independent variable (domain)
    
    failed = 0;
    % Interpolate and trucate NaN
    try [y,shift,endI] = interp_and_truncate_nan(Y);
    catch err
        failed = 1; shift = 1; endI = numel(Y);
    end
    
    consec_nan = find_consecutive_logical( isnan(Y(shift:endI)) );
    
    % Reject curves without the requisite number of non-NaN
    % data points. Also reject if too many consecutive NaN
    if max(consec_nan) < opt.nan_consec_thresh && any(Y > 0) ...
        && ~failed && numel(Y(~isnan(Y))) > opt.nan_thresh
        
        t = t( shift : shift + numel(y) - 1);
        
        % Establish the lower constraints
        lb = [0; t(1); opt.sigma_lb];
        % Establish the upper constraints
        ub = [nanmax(y); t(end) - opt.end_tol; opt.sigma_ub];
        
        % Perform multiple-fitting
        [gauss_p, residuals, J] = iterative_gaussian_fit( ...
            y,t,opt.alpha,lb,ub,opt.bg);
        
        % Turn on the 'fitted' flag
        this_cell.flag_fitted = 1;
        this_cell.opt = opt;
        
        % --- Get fitted curves, except for array of Gaussians ---
        %                     curve = synthesize_gaussians(gauss_p(:,2:end),t); % Get gaussians/time
        background = lsq_exponential(gauss_p(:,1),t); % Get background
        P = plot_peak_color(gauss_p(:,2:end),t); % Get colorized fit
        this_cell.fit_colorized = P; this_cell.fit_bg = background;
        this_cell.fit_time = t;
        this_cell.residuals = residuals;
        this_cell.jacobian = J;
        this_cell.params = gauss_p;
        this_cell.raw = y;
        
        % --- Get pulse info ---
        this_cell.num_fits = size(gauss_p,2) - 1;
        this_cell.fitID = [];
        
        fit_gausses = nan( numel(t), max(size(gauss_p,2) - 1, 1) );
        
        for j = 2 : size(gauss_p,2)
            
            % Assign fitID
            fitID = fitID + 1;
            
            embryo_fitID = fitID;
            
            fit = [fit Fitted(this_cell, gauss_p(:,j), embryo_fitID, opt)];
            
            % Construct Fitted object
            this_cell.fitID = [this_cell.fitID embryo_fitID];
            
            % Collect the cell-centric fitted curve for this peak
            fit_gausses(:,j - 1) = lsq_gauss1d(gauss_p(:,j),this_cell.fit_time);
        end
        
        this_cell.fit_gausses = fit_gausses;
    else % This cell will NOT be fitted
        this_cell.fit_colorized = NaN;
        this_cell.fit_bg = NaN;
        this_cell.fit_gausses = NaN;
        this_cell.fit_time = NaN;
        this_cell.residuals = NaN;
        this_cell.jacobian = NaN;
        this_cell.num_fits = NaN;
        this_cell.fitID = NaN;
        this_cell.num_tracks = NaN;
        this_cell.trackID = NaN;
        this_cell.params = NaN;
        this_cell.opt = NaN;
    end
    
%     new_cells(i) = this_cell;
    last_embryoID = this_cell.embryoID;
    
    display(['Done with cell #', num2str(stackID)])
    
end

end % fit_gaussians