function pulse = createFitFromTrack(pulse,trackID,opt)
%@Pulse.createFitFromTrack Using the stackID/embryoID and timing to create an
% artificial 'fit'.
%
% USAGE: pulse = pulse.createFitFromTrack(cells,trackID,fit_opt)
% xies@mit.edu Feb 2013

cells = pulse.cells;

% Extract track / make sure it's not duplicated
track = pulse.get_trackID(trackID);
c = pulse.get_cellID(track.cellID);
if isempty(track), display('Cannot create FIT: No track with trackID found.'); return; end

[cx,cy,ct] = pulse.get_xyt(track);

if isfield(pulse.changes,'fitsMadeFromTrack');
    foo = [pulse.changes.fitsMadeFromTrack.tracks];
    M = [[foo.cx];[foo.cy];[foo.cy]]';
    if any( bsxfun(@eq, M, [cx cy ct]) )
        display(['Track ID' num2str(trackID) ' already used to add a fit.']);
        return
    end
end

display(['Creating fit from trackID ' num2str(trackID)])

% Load already manually fitted params
try already_done = ...
        csvread( [ fileparts(pulse.tracks_mdf_file), '/', 'manual_fits.csv' ] );
catch err
    if strcmpi(err.identifier,'MATLAB:csvread:FileNotFound')
        already_done = [NaN NaN NaN];
    else
        rethrow(err);
    end
end

% check if already done
I = find( ...
    abs( already_done(:,1) - cx) < 2 ...
    & abs( already_done(:,2) - cy) < 2 ...
    & abs( already_done(:,3) - ct) < 10 ...
    );
if ~isempty(I)
    if numel(I) > 1
        keyboard;
    end
    % If already_done, then load recorded change
    params = already_done( I , 4:6 );
else
    % Launch the manual fit GUI
    params = manual_fit( ...
        [nanmean(track.dev_time) 20],c);
end

% Construct a new FITTED object from parameters
new_fit = Fitted( c, params, pulse.next_fitID, opt);
pulse.next_fitID = pulse.next_fitID + 1; %increment fitIDs

% Add into stack
[fits,errorflag] = pulse.fits.add_fit(new_fit);
if errorflag, return; end
fits(end).manually_added = 1;

pulse.fitsOI_ID = [pulse.fitsOI_ID fits(end).fitID];
pulse.fits = fits;

% Redo match/categorizing
pulse = pulse.match_pulse(pulse.match_thresh);
pulse = pulse.categorize_mapping;

% Update cell tracklist
% pulse.cells( [pulse.cells.stackID] == track.stackID ) = ...
%     pulse.cells( [pulse.cells.stackID] == track.stackID).addFit( pulse.fits(end));
c.addFit(pulse.fits(end));

% Record changes
%             this_change.fitID = pulse.fits(end).fitID;
%             this_change.trackID = trackID;
[cx,cy,ct] = pulse.get_xyt(track);
this_change.fits.cx = cx;
this_change.fits.cy = cy;
this_change.fits.ct = ct;
this_change.tracks.cx = cx;
this_change.tracks.cy = cy;
this_change.tracks.ct = ct;

if isfield(pulse.changes,'fitsMadeFromTrack')
    pulse.changes.fitsMadeFromTrack(end+1) = this_change;
else
    pulse.changes.fitsMadeFromTrack = this_change;
end

end % createFitFromTrack