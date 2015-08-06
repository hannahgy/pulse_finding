function varargout = movie(pulse, fitID, embryo_stack)
% MOVIE - Wrapper for MAKE_CELL_IMG to make a movie of a single
% fitted pulse (Fitted).
%
% USAGE: pulse.movie(fitID,embryo_stack);
% xies@mit.edu

% if nargin < 3
%     cells = embryo_stack;
%     embryo_stack = fitID;
%     fitID = fits.fitID;
% end

assert(numel(pulse) == 1)

% Extract this fit
this_fit = pulse.get_fitID(fitID);
if isempty(this_fit), varargout{1} = []; return; end
% Find its embryo
this_embryo = embryo_stack( this_fit.embryoID );
% The frames2load are the width-frames (1 sigma away)
frames = this_fit.margin_frames;
% Get vertices
h.vx = this_embryo.vertex_x;
h.vy = this_embryo.vertex_y;
% Get IDs
h.cellID = this_fit.cellID;
h.input = this_embryo.input;
% Get frames (no need to correct for t0, done in
% MAKE_CELL_IMG)
h.frames2load = frames;

h.channels = {'Membranes','Myosin'};

% Pad the curve
this_cell = pulse.get_cellID( this_fit.cellID );
h.measurement = nan(size( this_cell.dev_time ));
h.measurement( this_fit.margin_frames ) = this_fit.fit;

% Turn on segmentation border
h.border = 'off';

%             figure
F = make_cell_img(h);

if nargout, varargout{1} = F; end

end %movie
