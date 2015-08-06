function obj = find_pulse_by_xyt(pulse,obj_type,cx,cy,ct)
%FIND_PULSE_BY_XYT
% Use centroids, and mean dev_time to robustly find the
% specified object: 'track' or 'fit'

c = pulse.cells;
cframe = findnearest(ct,c(1).dev_time);
if numel(cframe) > 1, cframe = cframe(1); end
x = cat(2,c.centroid_x);
y = cat(2,c.centroid_y);

x = x(cframe,:); y = y(cframe,:);
d = (x-cx).^2 + (y-cy).^2;

[~,which] = min(d);
% cellID = c(which).cellID;

switch obj_type
    case 'track'
        
        obj = pulse.find_tracks_from_cell(c(which));
        which = findnearest(cellfun(@nanmean,{obj.dev_time}),ct);
        obj = obj(which);
        
    case 'fit'
        
        obj = pulse.find_fits_from_cell(c(which));
        which = findnearest(cellfun(@nanmean,{obj.dev_time}),ct);
        obj = obj(which);
        
    otherwise
        error('Can only take ''type'' or ''fit'' as object type')
end

if isempty(obj), obj = []; end

end