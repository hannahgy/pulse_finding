classdef Pulse
    %--- PULSE ------------------------------------------------------------
	%A housekeeping class to keep track of TRACK pulses and FITTED pulses,
	% as well as the bi-directional correspondence between the two sets.
	%
	% Matches tracked pulses (TRACK) to fitted pulses (FITTED), with a generalized
	% class MatchTrackFit to handle the two-way mapping between the two sets of
	% objects, which could contain the following categories:
	%	1) one-to-one matches
	%	2) missed objects from either set
	%	3) objects from one set mapped onto the same object in the other.
	%
	% Methods are implemented to handle manual adjustment to the correspondence
	% to resolve missed objects, as well as merged objects. To this end a graphical
	% output is also implemented under @pulse.graph.
	% 
	% Properties (private)
	%	fits - the complete array of all fits, does not have to exclude non-tracked
	%		cells/embryos
	%	fitsOI_ID - the fitID (unique ID code for each fitted pulse) of interest, i.e.
	%		those from cells with tracked pulses
	%	tracks - the array of TRACK pulses for this embryo
	%	tracks_mdf_file - the MDF filename
	%	cells - an array of CELLOBJ of data from cells found in this embryo
	%
	% Properties (public)
	% 	embryoID - the index of this embryo (see also LOAD_EDGE_SCRIPT)
	%	changes - a structure containing information about the manual changes
	%		performed on the PULSE object
	%
	% Methods
	% 
	% --- Construction methods ---
	%	Pulse - constructor object
	%	.match_pulse - performs matching between .tracks and .fits, and constructs
	%		the .map property from the resulting two-way mapping
	%	.categorize_mapping - creates/updates the .categories property from the .map
	%		includes the following fields:
	%		.one2one - one-to-one matches
	%		.merge - two Track per Fitted
	% 		.split - one Track per two Fitted
	%		.miss - Track with no Fitted
	% 		.add - Fitted with no Track
    %   .search_catID Search for the category_index (catID) of a pulse
    %       (track/fit).
    %   .cat - concatenate two Pulse objects (for example from different
    %       embryos)
    % 
    % --- toArray methods ---
    %   getCells - return all CellObjs from a Pulse vector as an array
    %   getFits - return all Fitted as array
    %   getTracks - return all Track as array
    %
    % --- Individual access ---
    %   get_cellID - access individual(s) by cellID
    %   get_fitID - access individual(s) by fitID
    %   get_trackID  - access individual(s) by trackID
    %   find_pulse_by_xyt - find Fitted/Track by spatiotemporal coordinate
    % 
    % --- Find fits/track from cell and vice versa ---
    %   find_cells_with_fit
    %   find_cells_with_track
    %   find_fits_from_cell
    %   find_tracks_from_cell
    %
    % --- Putting data into CellObj/Fitted ---
    %   assign_datafield - put data into a field in Fitted array
    %   align_fits - align Fitted by Gaussian center
    %   bin_fits - assign the intra-embryo percentile rank
    %   interpolate_traces - interpolate Fitted traces so they have the same
    %      temporal resolution
    %   retrace - retrace out subsequences with different temporal widths
    %      around Fit centers
    %   measure_fits - populate myosin and area data from CellObj into
    %      Fitted arrays
    %   get_cell_measurement - return a matrix of CellObj fields
    %
    % --- Making measurements on CellObj ---
    %   get_frequency - returns the frequency of pulsing
    %   estimate_pulsing_params - returns a gamma distribution fit of
    %      pulsing frequency and number of cell-neighbors count
    %   make_binary_sequence - returns a binary timeseries of whether a
    %      given cell is pulsing or not
    %   get_pulsing_trajectories - returns the temporal sequence of a
    %      cell's pulses and their respective behaviors
    %   get_pulse_transition_grpah - returns the nodes and edges of a
    %      bipartite graph for transition b/w pulse behaviors
    %   get_pulse_transition_matrix - same as above but in matrix form
    %
    % --- Making measurements on Fitted ---
    %   get_fit_measurement - return a matrix of Fitted fields (no other
    %      processing is done)
    %   get_corrected_measurement - like .get_fit_measurement but also does
    %      interpolation so all pulses have the same temporal resolution
    %   get_myosin_persistence - return myosin_persistence for all Fitted
    %   percent_overlap - measures the fraction of Fitted subsequences that
    %      overlap with each other
    %   fcm_cluster - perform FCM clustering
    %   fcm_stability - returns Rand Index metrics of FCM cluster stability
    %
    % --- Neighbor-neighbor (fit or cell) analysis ---
    %   find_near_fits - find fits in the spatiotemporal neighborhood
    %       (user-defined) as other fits
    %   find_non_edge - finds fits that are not in the border of
    %      segmentation
    %   get_num_near - returns the # of fits nearby other fits
    %   bootstrap_cluster_label - bootstrap test for cluster label
    %   simulate_pulsing - generates random pattern of pulsing
    %
    % --- Matching tracks to Fits ---
    %   match_tracks_to_fits - initial matching
    %   categorize_mapping - categorizes mapping as one2one, missing, etc.
    %   search_catID - returns map elements
    %
	% --- Manual editing methods ---
    %   rename_embryoID - consistently renames the embryoID
	%	search_catID - given a tracked or fit from a PULSE object, find the index
	%		of that object within its current .cagetory
	%	removePulse - remove a track/fit and update the mapping
	%	createFitFromTrack - create an artificial FITTED from roughly where a TRACK
	%		was... uses MANUAL_FIT
	%	createTrackFromFit - create an artificial TRACK from roughly where a FITTED
	%		was...
	%	reassignFit - re-assign a FITTED to a TRACK, will only work if neither have
	% 		prior assignments
	%	read_changes - given a .change structure, edit current Pulse object
    %
    % --- Edit embryo-level paramsters ---
    %   rename_embryoID - reindex all relevant IDs
    %   adjust_centers - adjust the reference time used to construct
    %       dev_time
    %
    % --- Saving methods ---
    %   export_manual_fits - writes manually fitted parameters into a CSV
    %       file1
	%	export_changes - writes all changes into a CSV
    %   get_xyt - returns the spatiotemporal coordinate of a Fitted/Track
    %
	% --- Display methods ---
	%	graph - Generates a 1x3 subplot of the TRACK/FIT/CELL
	% 	display - In-line display, reporting the number of objects and the quality of
	%		matching
    %
    % ------ Display Fitted ---
    %   plot_binned_fits - plot average myosin + area for each bin
    %   plot_heatmap - returns myosin and area heatmap for all Fitted
    %   plot_single_pulse - plot a single Fitted (myosin + area)
    %   movie - makes a movie of a single Fitted
    %
    % ------ CellObj display ---
    %   plot_cells_aligned - plot avg property (e.g. area) for all cells
    %      found in a pulse
    %
    % See also: CELLOBJ, FITTED, TRACK, FIND_ONE2ONE
	%
	% xies@mit.edu April 2013.

    properties
        
        fits 		% Complete set of FITTED pulses from ALL cells, including non-tracked and other embryos
        fitsOI_ID 	% fitID from only trakced cells in TRACKS
		fit_opt 	% The fitting option file for this embryo
        tracks 		% Set of TRACK pulses for this embryo
        tracks_mdf_file % The filename of the MDF file from which .tracks was loaded
        cells 		% An array of CELL of cells in this embryo (Contains the raw data)
        input      % input array with info about identifiy of each embryo
        next_fitID  % bookkeeping
        
        map			% The two-way mapping between Track and Fit
        match_thresh % The threshold of frames overlap above which a TRACK and a FITTED is matched (usually 1)
        categories %  A structure containing the different categories of matches
        
        embryoID
        changes
        bootstrapped
        
    end
    
    methods % Dynamic methods
        
% --------------------------- Constructor ---------------------------------

        function pulse = Pulse(fits,cells,opt,input)
			%PULSE Constructor for the Pulse object (see main documentation)
			% Will not generate the .map property.
			%
			% USAGE: pulse = Pulse(fits,cells,opts,input);
            
            if nargin > 0 % empty case (for object arrays)
                
                assert( unique([fits.embryoID]) == unique([cells.embryoID]) , ...
                    'CellObj and Fitted arrays should be from the same embryo.' );
%                 fitsOI_ID = [fits(ismember( [fits.cellID],[tracks.cellID] )).fitID];
                
                pulse.fits = fits.copy();
                pulse.cells = cells.copy();
                pulse.fit_opt = opt;
                pulse.next_fitID = fits(1).embryoID * 100000;
                pulse.input = input;
                pulse.embryoID = cells(1).embryoID;
            
            end % non-empty constructor
            
        end % Constructor
        
% ------------------- Cell/Fitted accessing -------------------------------

        % Gives access to all arrays
        function pulse = get_embryoID(pulse,embryoID)
            pulse = pulse(ismember([pulse.embryoID],embryoID));
        end
        
        % Convert to array
        function fits = getFits(pulse), fits = [pulse.fits]; end
        function cells = getCells(pulse), cells = [pulse.cells]; end
        function tracks = getTracks(pulse), tracks = [pulse.tracks]; end
        
        function fits = get_cluster(pulse,label)
            % Returns the cluster behavior
            % USAGE: filtered = fits.get_cluster(1:3)
            %   ABOVE will return all fits with cluster label 1-3.
            fits_array = [pulse.fits];
            fits = fits_array( ismember([ fits_array.cluster_label ], label) );
        end
        
        % Access individuals by ID
        fits = get_fitID(pulse,fitID);
        cells = get_cellID(pulse,cellID);
        tracks = get_trackID(pulse,trackID);
        % Find from fits/tracks to cells and vice versa
        fits = find_cells_with_fit(pulse,fits);
        cells = find_fits_from_cell(pulse,cells);
        tracks = find_cells_with_track(pulse,cells);
        cells = find_tracks_from_cell(pulse,tracks);
        first_fits = get_first_fit(pulse);
        
        % Access Fit/track by spatiotemporal coordinate
        obj = find_pulse_by_xyt(pulse,obj_type,cx,cy,ct);
        
% ------------------- Cell/Fitted data pre-processing ---------------------
        
        assign_datafield(pulse,data,name);
        align_fits(pulse,name,measurement);
        interpolate_traces(pulse,name);
        retrace(pulse, opts);
        measure_fits(pulse)
        bin_fits(pulse,range)
        A = get_cell_measurement(pulse,measurement_name);
        
% ---------------------- Cell-level analysis ------------------------------

        % Pulsing analysis
        [freq,center] = get_frequency(pulse);
        [freq,neighbor_count] = estimate_pulsing_params(pulse);
        binary = make_binary_sequence(pulse);
        [adj,nodes] = get_pulsing_trajectories(pulse);
        [adj,nodes] = get_pulse_transition_graph(pulse);
        W = get_pulse_transition_matrix(pulse);
        
%         [fits_bs,cells_bs] = monte_carlo_stackID(pulse)
        
% ------------------------ Pulse measurements -----------------------------
        
        A = get_fit_measurement(pulse,measurement_name);
        M = get_corrected_measurement(pulse,meas,input);
        myosin_persistence = get_myosin_persistence(fits);
        [perc,varargout] = percent_overlap(fits,cells);
        [avgRI,stdRI,avgR_random,stdRI_random] = fcm_stability(pulse,ks2try);
        A = fcm_cluster(pulse,k,datafield,max_nan);
        
% -------------------- Neighboring fits measurements ---------------------
        
        fits = find_near_fits(pulse,neighbor_def);
        f = find_non_edge(fits,cells);
        num_near = get_num_near(pulse,neighbor_definition,window);
        fits_bs = bootstrap_cluster_label(fits);
        [fits_bs,cells_bs] = simulate_pulsing(fits,cells,freqHat);

% --------------------------- Array handling ------------------------------
        
        pulse = horzcat(pulse1,pulse2);
        
% --------------------------- Mapping -------------------------------------
        
        pulse = match_tracks_to_fits(pulse,tracks,track_filename,threshold)
        pulse = categorize_mapping(pulse);
%         pulse = match_pulse(pulse,threshold);
        catID = search_catID(pulse,type,pulseID);
        
% --------------------- Edit pulse/tracks ---------------------------------
        
        pulse = removePulse(pulse,type,pulseID);
        pulse = createTrackFromFit(pulse,fitID);
        pulse = createFitFromTrack(pulse,trackID,opt);
        pulse = reassignFit(pulse,fitID,newTrackID);
		pulse = read_changes( pulse, changes );
        
% ----------------------- Saving / exporting ------------------------------

        export_manual_fits(pulse);
        [cx,cy,ct] = get_xyt(pulse,fit);
        
        function export_changes( pulse )
            %EXPORT_CHANGES
            % Export all .changes to a .mat file
            changes = pulse.changes;
            save( [fileparts(pulse.tracks_mdf_file), '/', 'changes.mat'], 'changes');
            % export manual changes
            pulse.export_manual_fits;
        end % export_changes
        
% ---------------------- graph/display ------------------------------------
        
        % Display for Pulse itself
        varargout = graph(pulse,cat,ID,axes_handle);
        disp(pulse);
        
        % Fitted display
        plot_binned_fits(fits);
        plot_heatmap(fits,sortname);
        fig = plot_single_pulse(fit,fitID);
        varargout = movie(pulse, fitID, embryo_stack)
        
        % CellObj alignment
        plot_cells_aligned(pulse,name2plot);
        
% ---------------------- Edit embryo-level parameters ---------------------

        function pulse = adjust_dev_time(pulse,input)
            % ADJUST_DEV_TIME
            % 
            % Adjust the .center and .dev_time of CellObj and Fitted
            % objects with respect to new reference development itme.
            
            if numel(pulse) > 1
                error('Only one embryo please.')
            end
            
            old_tref = find(pulse.cells(1).dev_time == 0);
            
            pulse.input = input;
            new_tref = input.tref;
            dt = input.dt;

            pulse.cells.adjust_dev_time(old_tref,new_tref,dt);
            pulse.fits.adjust_dev_time(old_tref,new_tref,dt);
            
        end % Adjust adjust_dev_time
        
        function pulse = rename_embryoID(pulse,embryoID)
            % Rename all the FITTED and CELLOBJ from an old embryoID into a
            % new embryoID.
            %
            % USAGE: pulse = pulse.rename_embryoID(newID)
            %
            %@TODO:  ALSO UPDATE the MAP!!!!!
            
            old_embryoID = pulse.embryoID;
            pulse.embryoID = embryoID;
            pulse.tracks = pulse.tracks.rename_embryoID(embryoID);
            pulse.fits = pulse.fits.rename_embryoID(embryoID);
            pulse.cells = pulse.cells.rename_embryoID(embryoID);
            
            for i = 1:numel(pulse.fitsOI_ID)
                
                fID = pulse.fitsOI_ID(i);
                base = 10.^floor(log10(fID) - log10(old_embryoID));
                fID = fID - old_embryoID*base + embryoID*base;
                pulse.fitsOI_ID(i) = fID;
                
            end
            
        end % rename_embryoID
        
        function pulse_bs = clear(pulse)
            pulse_bs = pulse;
            for e = 1:numel(pulse)
                pulse_bs(e).fits = pulse(e).fits.copy.clearCell;
                pulse_bs(e).cells = pulse(e).cells.copy.clearFitsTracks;
            end
        end
        
    end % Dynamic methods
    
end
