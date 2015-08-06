classdef CellObj < handle
	%--- CellObj ----------------------------------------------------------
    % Collects EDGE measurements for each cell, as well as the TRACK
	% and FITTED pulses found in that cell.
	%
	% PROPERTIES
    %   embryoID
	%   cellID
	%   stackID
    %   folder_name
    %   
    %   area
    %   area_sm
    %   centroid_x
    %   centroid_y
    %   identity_of_neighbors_all
    %   vertex_x
    %   vertex_y
    %   anisotropy_xy
    %   anisotropy
    %   myosin_intensity_fuzzy
    %   myosin_sm
    %
    %   dev_frame
    %   dev_time
    %
    % Properties (public)
    %   flag_fitted
    %   fit_colorized % RGB colorization of multiple peaks, for colored movie display (see MAKE_PULSE_MOVIE)
    %   fit_bg		% Background fit
    %   fit_gausses	% Sum of all gaussians fitted
    %   fit_time	% The time-domain of fitted pulses
    %   raw         % Raw curve
    %   residuals   % Residuals
    %   jacobian    % Jacobian
    %   params      % Keep track of pulse parameters
    %   num_fits    % Number of pulses found
    %   fitID       % The fitIDs of FITTED found in this cell
    %   opt         % fit_opts
    %
    % Methods
	%	--- Constructor ---
	%		CellObj - use from embryo2cell
	%	--- Fit pulses ---
	%		fit_gaussians - returns updated CellObj array
	%			and FITTED array
	%	--- Edit fit/tracks ---
	%		addFit - add fit to a cell
	%		removeFit - remove a fit from a cell
	%		addTrack - add a track to a cell
	%		removeTrack - remove trackID from a cell's records
    %       clearFitsTracks -- remove all fits and track from a cell's
    %                          records
	%	--- Array set/access ---
	%		get_stackID - search by stackID
	%		get_fitID - search for cells containing fitID
	%		get_trackID - search for cells containing trackID
    %       get_embryoID - only return embryos of given embryoID
	%		get_embryoID_cellID - search using embryoID + cellID (useful
    %           coming from EDGE)
    %       adjust_dev_time - adjust dev_time to reflect new reference time
    %       get_nearby - get cells nearby within a given radius
	%	--- Visualization/display ---
	%		make_mask - returns a binary BW image of the cell
	%		visualize - plots myosin + area v. time
	%		movie - makes movie of cell
    %   --- Analysis ---
    %       get_frequency
    %       get_pulsing_trajectories
    %       get_pulsing_transition_matrix
    %       get_pulsing_transition_graph
    %       get_adjacency_matrix
    %       get_neighbor_angle
    %       get_coronal_measurement
    %       monte_carlo_stackID
    %       make_binary_sequence
	%
	% See also: PULSE, FITTED, TRACK, EMBRYO2CELL
	%
	% xies@mit.edu April 2013.
    
    properties
        
        % IDs
        embryoID % embryo index
        cellID	% cellID in EDGE
        stackID	% stackID in collated cell stack
        folder_name % EDGE folder name
        
        % measurements
        area		% area time series
        area_sm		% smoothed area
        centroid_x	% x-centroid of cell
        centroid_y	% y-centroid of cell
        identity_of_neighbors_all	% cell arrays of neighboring stackIDs
        myosin_intensity % unfuzzy
        myosin_intensity_fuzzy		% myosin intensity (fuzzy border)
        myosin_sm					% smoothed fuzzy myosin intensity
        vertex_x	% vertices, x-coordinates
        vertex_y	% vertices, y-coordinates
        anisotropy  % anisotropy (shape)
        anisotropy_xy % anisotropy (projection)
        
        % Time
        dev_time	% multiple-embryo-aligned, developmental time
        
        % Fitting parameters / statistics
        flag_fitted	% Flagged if not skipped in fitting
        fit_colorized % RGB colorization of multiple peaks, for colored movie display (see MAKE_PULSE_MOVIE)
        fit_bg		% Background fit
        fit_gausses	% Sum of all gaussians fitted
        fit_time	% The time-domain of fitted pulses
        raw         % Raw curve
        residuals   % Residuals
        jacobian    % Jacobian
        params      % Keep track of pulse parameters
        num_fits    % Number of pulses found
        fitID       % The fitIDs of FITTED found in this cell
        opt         % fit_opts
        
        % Tracking
		flag_tracked % Flagged if tracked
        num_tracks  % number of tracks found in this cell
        trackID     % trackID of tracks found in cell
        
    end % Public properties - to do with track/fit
    
    methods
        
        function obj = CellObj(this_cell)
			% CellObj - Contructs object of CellObj class. Use from EMBRYO2CELL
            if nargin > 0
                field_names = fieldnames(this_cell);
                for i = 1:numel(field_names)
                    obj.(field_names{i}) = this_cell.(field_names{i});
                end
                obj.flag_tracked = 0;
                obj.num_tracks = 0;
            end
        end % Constructor
        
        function b = copy(a)
            %COPY Make a deep copy of a Fitted object (array supported).
            num_cells = numel(a);
            if num_cells == 0
                b = Fitted;
                return
            end
            
            b(1,num_cells) = CellObj; % Do not use 1:num_fits to initiate
            props = properties( a(1) );
            for i = 1:num_cells
                for p = ensure_row(props)
                    b(i).(p{:}) = a(i).(p{:});
                end
            end
        end
        
        
% ---------------------- Editing fit/tracks -------------------------------
        
        fit = fit_gaussians(cells,opts);
        
        cellobj = addFit(cellobj,fit);
        cellobj = removeFit(cellobj,fitID);
        cellobj = addTrack(cellobj,trackID);
        cellobj = removeTrack(cellobj,trackID);
        cell_array = clearFitsTracks(cell_array)
        
% --------------------- Array set/access ----------------------------------
        
%         function obj = get_stackID(obj_array, stackID)
%             %@Cell.get_stackID Returns the obj from an array with the given
%             % stackID
%             obj = obj_array( ismember([ obj_array.stackID ],stackID ));
%         end % get_stackID
        
%         function obj = get_fitID(obj_array, fitID)
%             %@Cell.get_fitID Returns the obj from an array with the given
%             % fitID
%             obj = obj_array(...
%                 cellfun(@(x) (any(x == fitID)),{obj_array.fitID}) );
%         end % get_fitID

        fits = getFits(cells,fits);
        
        function obj = get_trackID(obj_array, trackID)
            %@Cell.get_fitID Returns the obj from an array with the given
            % trackID
            obj = obj_array([obj_array.trackID] == trackID);
        end % get_trackID
        
%         function obj = get_embryoID(obj_array,embryoID)
%             obj = obj_array(ismember([obj_array.embryoID],embryoID));
%         end % get_embryoID
        
%         function obj = get_embryoID_cellID(obj_array,embryoID,cellID)
%             obj = obj_array( ...
%                 [obj_array.embryoID] == embryoID & ...
%                 ismember([obj_array.cellID], cellID) ...
%                 );
%         end % get_embryoID_cellID
        
        function obj = get_curated(obj_array)
            obj = obj_array([obj_array.flag_tracked] == 1 & ...
                [obj_array.flag_fitted] == 1);
        end
        
        nearby_cells = get_nearby(obj_array,stackID,radius,reference_frame);
        update_measurements(cells,embryo_stack);
        
        function adjust_dev_time(cells, old_tref, new_tref, dt)
            %ADJUST_DEV_TIME Recalculate dev_time according to new .tref
            % Properties that will be adjusted:
            %   dev_time
            
            % Only one embryo
            if numel(unique([cells.embryoID])) > 1
                error('Cells from only one embryo please.');
            end
            
            for i = 1:numel(cells)
                this_cell = cells(i);
                this_cell.fit_time = this_cell.fit_time - ...
                    (new_tref - old_tref)*dt;
                this_cell.dev_time = this_cell.dev_time - ...
                    (new_tref - old_tref)*dt;
                cells(i) = this_cell;
            end
        end % adjust_dev_time
        
        function cells = rename_embryoID(cells,embryoID)
            % Rename all cells into a new embryoID
            % Please use from PULSE only to ensure FIT objects are
            % similarly renamed
            
            old_embryoID = cells(1).embryoID;
            [cells.embryoID] = deal(embryoID);
            
            for i = 1:numel(cells)
                fIDs = cells(i).fitID;
                tIDs = cells(i).trackID;
                % regular fits are base 1000 while manual ones are base
                % 10000
                base = 10.^floor(log10(fIDs) - log10(old_embryoID));
                fIDs = fIDs - old_embryoID*base + embryoID*base;
                base = 10.^floor(log10(tIDs) - log10(old_embryoID));
                tIDs = tIDs - old_embryoID*(base) + embryoID*base;
                
                cells(i).fitID = fIDs;
                cells(i).trackID = tIDs;
                
                cells(i).stackID = 1000*embryoID + cells(i).cellID;
                
            end
        end % rename_embryoID
        
%---------------------- Visualization/display -----------------------------
        
        mask = make_mask(obj_array, frames, input);
        [x,y] = make_polygon(obj_array,t,input,filename);
        visualize(cell,handle);
        H = plot_aligned(cells,name2plot,varargin)
        M = movie(cells,stackID,embryo_stack);
        
% ------------------------- Analysis --------------------------------------
        
        % Tissue analysis
        N = get_adjacency_matrix(cells,method);
        angles = get_neighbor_angle(cellx,celly,frame);
        corona_measurement = get_corona_measurement(cells,measurement);
        % Pulsing analysis
        [fits,cells] = monte_carlo_stackID(cells,fits,opt);
        
   end % End methods 
   
end
