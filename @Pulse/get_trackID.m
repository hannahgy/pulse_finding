function tracks = get_trackID(pulse,trackID)
%GET_CELLID Return the Fitted(s) with the given fitID(s) from a single pulse.
% 
% USAGE: cells = pulse.get_trackID(trackID);

assert(numel(pulse) == 1)

tracks = pulse.tracks;
tracks = tracks( ismember([tracks.trackID],trackID) );

end