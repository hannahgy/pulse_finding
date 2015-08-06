function obj_array = clearFitsTracks(obj_array)
% Clear all records of fits and tracks from a cell
for i = 1:numel(obj_array)
    
    cellOI = obj_array(i).copy;
    % reset all attributes relating to track/fit
    cellOI.flag_fitted = 0;
    cellOI.flag_tracked = 0;
    cellOI.fit_colorized = [];
    cellOI.fit_bg = [];
    cellOI.fit_gausses = [];
    cellOI.fit_time = [];
    cellOI.raw = [];
    cellOI.residuals = [];
    cellOI.jacobian = [];
    cellOI.params = [];
    
    cellOI.num_fits = 0;
    cellOI.num_tracks = 0;
    
    cellOI.fitID = NaN;
    cellOI.trackID = NaN;
    
    obj_array(i) = cellOI;
    
end
end