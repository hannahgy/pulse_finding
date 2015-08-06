function tracks = load_mdf_track(mdf_matrix,embryo_struct,min_frame,cells)
%LOAD_MDF_TRACKS Loads tracked pulses from a MDF matrix into a TRACK
% object
%
% USAGE: tracks = load_mdf_track(mdf_matrix,embryo_struct,min_frame,cells)
%
% INPUT: MDF_MATRIX - output of READ_MDF
%        EMBYRO_STRUCT - 1 single structure of the EDGE output (see
%           LOAD_EDGE_DATA)
%        MIN_FRAME_OVERLAP - mininum overlap b/w track and fit to call a
%           match between them (e.g. 2 frame)
%        cells - CellObj array to map everything onto.
%
% See also TRACK, INPOLYGON
% xies@mit.edu

assert( numel(embryo_struct) == 1, 'Only 1 embryo_struct please.');
assert( embryo_struct.input.embryoID == unique([cells.embryoID]), ...
    'EmbryoID must be the same from embryo_struct to CellObjs.');

num_cells = numel(cells);

% Load time, cell number, and input from embryo_struct
dev_time = embryo_struct.dev_time;
% dev_frame = embryo_struct.dev_frame;
input = embryo_struct.input;

% Get rid of other zslices and squeeze out singleton
vx = squeeze( embryo_struct.vertex_x );
vy = squeeze( embryo_struct.vertex_y );
cx = squeeze( embryo_struct.centroid_x ) / input.um_per_px;
cy = squeeze( embryo_struct.centroid_y ) / input.um_per_px;

% Get the number of tracks in the matrix
num_tracks = max(unique(mdf_matrix(:,1)));
% Construct padded num_cell count for setting correct stackID
% num_cell_pad = cumsum([0 num_cells]);

% clear cells records
for i = 1:numel(cells)
    cells(i).trackID = [];
    cells(i).num_tracks = 0;
end

trackID = 0;
for i = 1:num_tracks

	% Extract one track from matrix
	this_mdf = mdf_matrix(mdf_matrix(:,1) == i,:);
	
	% If nothing found, move onto next number
	if isempty(this_mdf), continue; end

	% Get the active frames
	frames = this_mdf(:,end);
    frames = frames - input.t0;

	% Make sure there are more than a threshold number of active frames
	if numel(frames) < min_frame; end
	
	% Make sure that the frame doesn't lie beyond the well-tracked region
	if mean(frames) > input.last_segmented, continue; end
    % Make sure that the frame doesn't begin before EDGE
	if mean(frames) <= 0, continue; end
    
    frames( frames <= 0 ) = [];
    
	% Get the coordinates of the track and the centroid of track
	x = this_mdf(:,3);
	track_cx = mean(x);
	y = this_mdf(:,4);
	track_cy = mean(y);
    
    % Convert image-frames into dev-frames
%     f0 = find( 1 == dev_frame );
%     frames = frames + f0 - 1;

	% Sort the distance between track cenroid to all cell centroids
	dists = sqrt( ...
		(track_cx - cx(frames(1),:)).^2 + (track_cy - cy(frames(1),:)).^2);
	[~,order] = sort(dists);

	found = 0; index = 0;
	% Try until an overlapping cell is found using INPOLYGON
	while ~found && index < num_cells
		% Check if track center is within cell
		index = index + 1;
		found = inpolygon( track_cx,track_cy,...
			vx{frames(1),order(index)} , vy{frames(1),order(index)} );
	end
	if ~found, continue; end % Not in EDGEd cells

	% Assign trackID
	trackID = trackID + 1;
    this_track.trackID = trackID; % assign non-ambiguous trackID
    
	% Collect relevant information into track struct
    this_track.embryoID = input.embryoID; this_track.mdfID = i;
    this_track.cellID = order(index);
%     this_track.stackID = 1000*(embryoID) + order(index);

	% Collect the time/frame of track WRT aligned developmental time
    frames( frames > numel(dev_time) ) = [];
	this_track.dev_frame = ensure_row(frames);
	this_track.dev_time = dev_time(frames);
    
    % Update track info in CellObj (by reference)
	cells(order(index)).flag_tracked = 1;
    cells(order(index)).num_tracks = cells(order(index)).num_tracks + 1;
	% Construct Track object
    
	tracks(trackID) = Track(this_track);

end

% Filter out Tracks mapped onto non-fitted cells
unfit_cells = cells(~[cells.flag_fitted]);
tracks(ismember([tracks.cellID],[unfit_cells.cellID])) = [];

% for i = 1:numel(tracks)
%     tracks(i).trackID = i + input.embryoID*1000;
% %     cells( [cells.stackID] == tracks(i).stackID ).trackID = ...
% %         [cells( [cells.stackID] == tracks(i).stackID ).trackID i + input.embryoID*1000];
% end

end
