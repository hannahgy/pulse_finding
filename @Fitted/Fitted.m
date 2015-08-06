classdef Fitted < handle
    %--- FITTEd -----------------------------------------------------------
    % A fitted pulse as found by multiple-Gaussian fitting.
    %
    % Properties:
    %    (private)
    %       embryoID - embryo identifier
	%		cellID - EDGE ID from that embryo
	%		stackID - unique cell identifier across cells
	%		fitID - unique identifier for the FITTED object
	%		
	%		amplitude - height of pulse
	%		center - develpmental time corresponding to the maximum of the pulse
	%		width - sigma of the pulse
	%
    %       center_frame - the movie frame corresponding to the Gaussian
    %                  center
	%		margin_frames - standard margin of sub-frames set by fit_opt
	%		width_frames - sub-frames set by 1 sigma before and after center
	%		dev_time - developmental time corresponding to width_frames
	%		raw - the raw trace
	%		fit - the fitted trace
	%		residual
	%		aligned_time - time aligned across an array of pulses by their
	%			centers
	%		aligned_time_padded - NaN padded to make a matrix
	%		fit_padded - fitted trace padded to make a matrix
	%
	%		corrected_time - corrected for frame-rate differences
	%		myosin - myosin intensity
	%		myosin_rate - rate of myosin change
	%		area - cell area
	%		area_rate - rate of area change
	%		area_norm - cell area mean-subtracted
	%		anisotropy - anisotropy of cell shape
	%		corrected_myosin
	%		corrected_myosin_rate
	%		corrected_area
	%		corrected_area_rate
	%		corrected_area_norm
	%		measurement - placeholder for un-specified measurements
	%		corrected_measurement - placeholder
	%
	%	(public)
	%
	%		category - for use by PULSE
	%		bin - intra-embryo strength
	%		manually_added
	%		bootstrapped - flag if random shuffling was performed
	%	
	%		time_windows - neighborhood time-windows
	%		nearIDs - fitIDs of nearby fits, ordered WRT time_windows
    %		near_angles - angles b/w the fit and its neighbor listed in
    %   		nearIDs
	%		nearest_neighbor - fitID of nearest fit
    %       neighbor_cells
	%		
	%		cluster_label - should be in order
	%		cluster_weights - FCM weights
	%
    % Methods:
    %
    % --- Constructor ---
    %	Fitted - INPUT: CellObj, parameter, fitID, and opt_struct
    % --- Edit fit array ---
    %	add_fit - tries to append a given new_fit onto an array. Fails if
    %   	the same pulse already exists in the array
    %	removeFit - given a fitID, remove it from the object stack
    %   cat - concatenate two Fitted arrays
    % --- Comparator ---
    %   eq - right now will be equal if overlap of width_frame is > 3
    % --- Array access/set ---
    %   get_stackID - get all FITs with a certain stackID
    %   get_fitID - get the fit with the given fitID
    %   get_embryoID - get fit with the given embryoID
    %   clearCell - clear all records of which pulse belonged to which cell
    % --- Alignment ---
    %   adjust_dev_time - adjust Gaussian centers to tref
    %   get_corrected_measurement - temporarily puts given measurement and
    %      returns the aligned, interpolated .corrected_measurement
    % --- Array operations ---
    %   sort - sort according to field (default: amplitude)
    %   bin_fits - bin each embryo by amplitude
    % --- Visualization ---
    %   plot_binned_fits
    %   plot_heatmap (sorted)
    %   movie
    % --- Export ---
    %   export_field2csv (export given measurement to csv file)
    %   export_xyt - exports XYT of a set of pulses
    %
    % See also: PULSE, TRACK, CELLOBJ
    %
    % xies@mit.edu April 2013.
    
    properties
        
        % Initialized with
        % Identifying information
        embryoID
        cellID
%         stackID
        fitID
        
        % Pulse-specific data
        amplitude
        center
        width
        center_frame
        margin_frames
        width_frames
        dev_time
        raw
        fit
        residual
        aligned_time
        aligned_time_padded
        fit_padded
        opt % fitting options used
        
        % Added later
        corrected_time
        myosin
        myosin_rate
        area
        area_rate
        area_norm
        anisotropy
        corrected_myosin
        corrected_myosin_rate
        corrected_area
        corrected_area_rate
        corrected_area_norm
        % placeholder for further measurements
        measurement
        corrected_measurement
        
    end
    properties (SetAccess = public)
        
        category % one2one, added, missed, so forth
        bin % intra-embryo pulse strength bin
        % Flags to keep track during analysis
        manually_added % whether it's manually added
        bootstrapped % whether its randomly shuffled
        
        time_windows %
        nearIDs % fitIDs of 'nearby' fits, ordered w.r.t. time_windows
        near_angles % angle b/w the two pulses SEE: get_neighbor_angle
		nearest_neighbor % fitID of the nearest fit
        neighbor_cells % number of neighboring cells
        
        cluster_label
        cluster_weight
        
    end
    methods % Dynamic methods
        
% --------------------- Constructor ---------------------------------------
        
        function this_fit = Fitted(cell,params,fitID,opt)
            %Fitted Constructor - use from FIT_GAUSSIANS (array constructor)
            % Will populate the pulse-centric fields, like the margins, and
            % aligned curves and time.
            %
            % USAGE: fit = FITTED(cell,params,fitID)
            
            if nargin > 0
                % FitID
                this_fit.fitID = fitID;
                
                % Collect the relevant indices
                this_fit.embryoID = cell.embryoID;
                this_fit.cellID = cell.cellID;
%                 this_fit.stackID = cell.stackID;
                
                % Collect the parameters
                this_fit.amplitude = params(1);
                this_fit.center = params(2); this_fit.width = params(3);
                
                dev_time = cell.dev_time;
                num_frames = numel(dev_time);
                center_frame = findnearest(this_fit.center,dev_time);
                % ensure center_frame isn't two frames
                if numel(center_frame) > 1
                    center_frame = center_frame(1);
                end
                % Store center_frame
                this_fit.center_frame = center_frame;
                
                % Get pulse margin-time frame
                [left_margin,pad_l] = max([ center_frame - opt.left_margin , 1]);
                [right_margin,pad_r] = min([ center_frame + opt.right_margin , num_frames]);
                this_fit.margin_frames = left_margin:right_margin;
                
                % Get pulse width-time frame
                left_width = max( ...
                    center_frame - findnearest(this_fit.width,cumsum(diff(dev_time))) , 1);
                right_width = min( ...
                    center_frame + findnearest(this_fit.width,cumsum(diff(dev_time))) , num_frames);
                this_fit.width_frames = left_width:right_width;
                
                % Get pulse time WRT dev_time
                this_fit.dev_time = dev_time(left_width:right_width);
                
                % Collect the pulse-centric fitted curves
                x = dev_time( left_margin : right_margin );
                fitted_y = lsq_gauss1d( params , x );
                this_fit.fit = fitted_y;
                this_fit.raw = ...
                    ensure_row(cell.(opt.to_fit)( left_margin : right_margin ));
                res = fitted_y - this_fit.raw;
                this_fit.aligned_time = x - this_fit.center;
                
                % PAD the margin-time frame for plotting purposes
                if pad_l > 1
                    fitted_y = [ensure_row(nan(1 - (center_frame - opt.left_margin), 1)), fitted_y];
                    x = [ensure_row(nan(1 - (center_frame - opt.left_margin), 1)), x];
                    res = [ensure_row(nan(1 - (center_frame - opt.left_margin), 1)), res];
                end
                if pad_r > 1
                    fitted_y = [fitted_y, nan(1, (center_frame + opt.right_margin) - num_frames)];
                    x = [x, nan(1, (center_frame + opt.right_margin) - num_frames)];
                    res = [res, nan(1, (center_frame + opt.right_margin) - num_frames)];
                end
                this_fit.aligned_time_padded = x;
                this_fit.fit_padded = fitted_y;
                this_fit.residual = res;
                
                this_fit.manually_added = 0;
                this_fit.bootstrapped = 0;
                this_fit.opt = opt;
            end % From scratch
            
        end % constructor
        
        function b = copy(a)
            %COPY Make a deep copy of a Fitted object (array supported).
            num_fits = numel(a);
            if num_fits == 0
                b = Fitted;
                return
            end
            
            b(1,num_fits) = Fitted; % Do not use 1:num_fits to initiate
            props = properties( a(1) );
            for i = 1:num_fits
                for p = ensure_row(props)
                    b(i).(p{:}) = a(i).(p{:});
                end
            end
        end
        
% --------------------- Edit fit array -----------------------------------
        
        [obj_array,errorflag] = add_fit(obj_array,new_fit);
        obj_array = removeFit(obj_array,fitID);
        fits = rename_embryoID(fits,embryoID);
        
% --------------------- Comparator ----------------------------------------

        equality = eq(fit1,fit2)
        
% --------------------- Array access for Fitted/CellObj -------------------

        cells = findCell(fit,cells);
        
%         function fits = get_cellID(fit_array,cellID)
%             % Find the FIT(s) with the given cellID(s)
%             % usage: fitsOI = fits.get_cellID(cellID)
%             fits = fit_array( ismember([ fit_array.cellID ], cellID) );
%         end % get_stackID

%         function fits = get_embryoID(fit_array,embryoID)
%             % Find the FIT(s) with the given embryoID(s)
%             % USAGE: fitsOI = fits.get_embryoID(embryoID)
%             fits = fit_array( ismember([ fit_array.embryoID ], embryoID) );
%         end % get_embryoID
        
        fits = get_fitID(fit_array,fitID);
        
        fit_array = clearCell(fit_array);
        
% --------------------- Alignment functions -------------------------------
        
%         align_fits(fits,cells,name,measurement);
%         assign_datafield(fits,data,name);
%         interpolate_traces(fits,name,dt);
%         retrace(fits, cells, opts);
        adjust_dev_time(fits, old_tref, new_tref, dt);
        
%         M = get_corrected_measurement(fits,c,meas,input);
%         [cx,cy,ct] = get_xyt(fit,cell);
        
% --------------------- Array operations ----------------------------------
        
%         bin_fits(fits,range);
        fits_new = sort(fits,field);
        
% --------------------- Analysis ------------------------------------------
        
%         X = fcm_cluster(fits,k,datafield,max_nan);
%         
%         fits = find_near_fits(fits,cells,neighbor_def);
%         num_near = get_num_near(fits,cells,neighbor_definition,window);
%         
%         fits_bs = bootstrap_cluster_label(fits);
%         [fits_bs,cells_bs] = simulate_pulsing(fits,cells,freqHat);
%         
%         myosin_persistence = get_myosin_persistence(fits);
%         
% 	    [perc,varargout] = percent_overlap(fits,cells);
%         f = find_non_edge(fits,cells);
        
% --------------------- Visualization -------------------------------------

        varargout = movie(fits, fitID, embryo_stack, cells);
        
% ------------------------- Export ----------------------------------------
        
        export_field2csv(fits,filepath,fieldname);
%         [cx,cy,ct] = export_xyt( fits, cells, filename, traceback);
        
    end % Dynamic methods
    
end
