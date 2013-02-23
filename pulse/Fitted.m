classdef Fitted
    properties (SetAccess = private)
        
        % Initialized with
        embryoID
        cellID
        stackID
        fitID
        amplitude
        center
        width
        margin_frames
        width_frames
        img_frames
        dev_time
        raw
        fit
        aligned_time
        aligned_time_padded
        fit_padded
        
        % Added later
        myosin
        myosin_rate
        area
        area_rate
        area_norm
        corrected_myosin
        corrected_myosin_rate
        corrected_area
        corrected_area_rate
        corrected_area_norm
        
    end
    properties (SetAccess = public)
        
        category
        manually_added
        
    end
    methods % Dynamic methods
        function obj = Fitted(this_fit)
            % Constructor - use from FIT_GAUSSIANS (array constructor)
            if nargin > 0
                names = fieldnames(this_fit);
                for i = 1:numel(names)
                    [obj.(names{i})] = deal(this_fit.(names{i}));
                end
                obj.manually_added = 0;
            end
        end % constructor
        function obj_array = removeFit(obj_array,fitID)
            obj_array(obj_array.get_fitID(fitID)) = [];
        end
        
        function objs = get_stackID(obj_array,stackID)
            % Find the FIT(s) with the given stackID(s)
            objs = obj_array( ismember([ obj_array.stackID ], stackID) );
        end %get_stackID
        
        function objs = get_fitID(obj_array,fitID)
            % Find the FIT(s) with the given fitID(s)
            objs = obj_array( ismember([ obj_array.fitID ], fitID) );
        end %get_fitID
        
        function [fits] = align_fits(fits,measurement,name,opt)
            %ALIGN_PEAKS Aligns the global maxima of a given array of
            %FITTED objects
            % Will return also a given measurement aligned according to the
            % maxima. Updates the FITTED structure.
            %
            % SYNOPSIS: [fits,time] = align_peaks(fitss,cells,opt);
            
            num_fits = numel(fits);
            duration = numel(fits(1).margin_frames);
            l = opt.left_margin; r = opt.right_margin;
            
            center_idx = l + 1;
            
            for i = 1:num_fits
                
                frames = fits(i).margin_frames;
                %                 nonan_idx = fits(i).aligned_time_padded;
                %                 nonan_idx = ~isnan( nonan_idx );
                [~,max_idx] = max( fits(i).fit );
                
                left_len = max_idx - 1;
                
                m = nan(1, l + r + 1);
                m( (center_idx - left_len) : (center_idx - left_len + duration - 1) ) = ...
                    measurement( frames, fits(i).stackID );

                fits(i).(name) = m;
                
            end
            
        end % align_fits
        
        function [fits] = assign_datafield(fits,data,name)
            if size(data,1) ~= numel(fits)
                error('Data size must be the same as the number of FITTED objects.');
            end
            for i = 1:numel(fits)
                fits(i).(name) = data(i,:);
            end
        end % assign_datafield
        
        function pulses = resample_traces(pulses,name,dt,opt)
            %RESAMPLE_TRACES Uses INTERP1 to resample short traces
            %
            % [aligned_traces,aligned_time] = resample_traces(traces,embryoID,dt);
            %
            % xies@mit.edu Oct 2012
            
            traces = cat(1,pulses.(name));
            
            num_traces = size(traces,1);
            
            embryoIDs = [pulses.embryoID];
            if numel(embryoIDs) ~= num_traces
                error('The number of traces and the number of embryoID must be the same.');
            end
            
            aligned_dt = round(mean(dt)*100)/100;
            l = opt.left_margin; r = opt.right_margin;
            % w = floor(T/2);
            
            % aligned_traces = zeros([num_traces, l + r - 3]);
            aligned_t = (- l : r )*aligned_dt;
            
            % Resample using the SIGNAL_PROCESSING TOOLBOX
            for i = 1:num_traces
                pulses(i).(['corrected_' name]) = ...
                    interp1((-l:r)*dt(embryoIDs(i)),traces(i,:),(-(l-2):r-2)*aligned_dt);
                pulses(i).corrected_time = aligned_t;
            end
            
        end % resample_traces
        
    end
    
    methods (Static)
        
    end
    
end
