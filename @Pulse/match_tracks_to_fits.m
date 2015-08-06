function pulse = match_tracks_to_fits(pulse,tracks,track_filename,threshold)
%MATCH_TRACKS_TO_FITS Match FITTED objects to TRACK objects by generating a
% MatchTrackFit object, which is a wrapper for interacting with two
% container.Map objects (hashmaps).
%
% USAGE: pulse = match_tracks_to_fits(pulse,tracks,track_filename,threshold);

% Make sure that tracks are from the same embryo
assert(numel(pulse) == 1);
assert( unique([tracks.embryoID]) == pulse.embryoID, ...
    'Error making Pulse object: Tracks need to come from the same embryo as Pulse.');

% if strcmpi(class(tracks),'Track')
%     pulse.fits = fits;
%     pulse.tracks = tracks;
% else
%     pulse.fits = tracks;
%     pulse.tracks = fits;
% end

pulse.tracks = tracks;
fits = pulse.fits;

fitsOI_ID = [fits( ... filter out non-tracked cells
    ismember( [fits.cellID], [tracks.cellID] )).fitID];

pulse.tracks_mdf_file = track_filename;
pulse.fitsOI_ID = fitsOI_ID;

% save([fileparts(pulse.tracks_mdf_file) '/pulse_raw.mat'], 'pulse');

if isempty(pulse.map) || ...
        (~isempty(pulse.match_thresh) && pulse.match_thresh ~= threshold)
    pulse.categories = [];
    pulse.changes = [];
end

nbm = MatchTrackFit(pulse.tracks,pulse.fits,threshold);
pulse.map = nbm;
pulse.match_thresh = threshold;

end %
