function [cx,cy,ct] = get_xyt(pulse,fit)
%GET_XYT Return the x- and y-centroid location and timing of a pulse.
%
% USAGE: [cx,cy,ct] = pulse.get_xyt(fit);
%        [cx,cy,ct] = pulse.get_xyt(track);

switch class(fit)
    case 'Fitted'
        
        cell = pulse.find_cells_with_fit(fit);
        validateattributes(fit,{'Fitted'},{'scalar'});
        validateattributes(cell,{'CellObj'},{'scalar'});
        
        cframe = findnearest(nanmean(fit.dev_time),cell.dev_time);
        if numel(cframe) > 1, cframe = cframe(1); end
        
        cx = repnan(cell.centroid_x);
        cx = cx(cframe);
        
        cy = repnan(cell.centroid_y);
        cy = cy(cframe);
        
        ct = mean(fit.dev_time);
        
    case 'Track'
        
        track = fit;
        cell = pulse.find_cells_with_track(track);
        validateattributes(track,{'Track'},{'scalar'});
        validateattributes(cell,{'CellObj'},{'scalar'});
        
        cframe = findnearest(nanmean(track.dev_time),cell.dev_time);
        if numel(cframe) > 1, cframe = cframe(1); end
        cx = cell.centroid_x(cframe);
        cy = cell.centroid_y(cframe);
        if any(isnan([cx cy]))
            cframe = find_nearest_nonan(cell.centroid_x,cframe);
            
            cx = cell.centroid_x(cframe);
            cy = cell.centroid_y(cframe);
        end
        
        ct = nanmean(track.dev_time);
        
    otherwise
        error('Unrecognized input...')
        
end

end