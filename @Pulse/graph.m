function varargout = graph(pulse,cat,ID,axes_handle)
% Graph the selected cateogry
% USAGE: pulse.graph(category,ID,handles)
%        pulse.graph(category,ID)
%
% INPUT: category - string corresponding to category name, e.g.
%               'one2one' or 'merge'
%        ID - out of this category, a vector of IDs
%        axes_handle - subplot axes

% get the data
% fits = pulse.fits; tracks = pulse.tracks; cells = pulse.cells;

% find number of things to graph
category = pulse.categories.(cat);
num_disp = numel(ID);

% Default axes = gca
if nargin < 4, axes_handle = gcf; end

for i = 1:num_disp
    
    % obtain relevant highlight IDs (could be empty)
    fitID = category(ID(i)).fitID;
    trackID = category(ID(i)).trackID;
    % Get stackID
    if ~isempty(trackID)
        cellID = pulse.get_trackID(trackID(1)).cellID;
    else
        cellID = pulse.get_fitID(fitID(1)).cellID;
    end
    
    % Get time (for graphing
    cellOI = pulse.get_cellID(cellID);
    dev_time = cellOI.dev_time;
    dev_time = dev_time(~isnan(dev_time));
    
    % Extract fit/track of interest
    track = pulse.find_tracks_from_cell(cellOI); num_track = numel(track);
    fit = pulse.find_fits_from_cell(cellOI); num_fit = numel(fit);
    
    % --- Plot tracked pulses ---
    % handle subplots, plot to alternative parent if applicable
    h(1) = subplot(3, num_disp, i, 'Parent', axes_handle);
    
    %                 axes(h(1));
    binary_trace = concatenate_pulse(track,dev_time); % get binary track
    if ~isempty(trackID) % highlight pulse if applicable
        on = highlight_track(track,trackID);
        binary_trace(on,:) = binary_trace(on,:) + 3;
    end
    if num_track > 1 % Plot
        %                     if nargin > 4
        imagesc(dev_time,1:num_track,binary_trace,'Parent',h(1));
        
    elseif num_track == 1
        plot(h(1),dev_time,binary_trace);
    else
        cla(h(1));
    end
    set(h(1),'Xlim',[min(dev_time) max(dev_time)]);
    xlabel(h(1),'Develop. time (sec)');
    %                 title(h(1),['Manual: #' num2str(track(1).trackID)])
    
    % --- Plot fitted pulses ---
    
    h(2) = subplot(3, num_disp, num_disp + i, 'Parent', axes_handle);
    
    %                 axes(h(2));
    binary_trace = concatenate_pulse(fit,dev_time); % get binary track
    if ~isempty(fitID) % highlight pulse if applicable
        on = highlight_track(fit,fitID);
        binary_trace(on,:) = binary_trace(on,:) + 3;
    end
    if num_fit > 1 % Plot
        imagesc(dev_time,1:num_fit,binary_trace,'Parent',h(2));
    elseif num_fit == 1
        plot(h(2),dev_time,binary_trace);
    else
        cla(h(2));
    end
    set(h(2),'Xlim',[min(dev_time) max(dev_time)]);
    xlabel(h(2),'Develop. time (sec)');
    title(h(2),['Fitted'])
    
    % --- Plot cell raw data ---
    h(3) = subplot(3, num_disp, 2*num_disp + i, 'Parent', axes_handle);
    
    cellOI.visualize( h(3) );
    linkaxes( h , 'x');
    
end % End of for-loop

if nargout > 0
    varargout{1} = [pulse.find_tracks_from_cell(cellOI).trackID];
    varargout{2} = [pulse.find_fits_from_cell(cellOI).fitID];
end

% --- Sub functions ---
    function binary = concatenate_pulse(pulse,time)
        % Creates a binary track of all pulses
        binary = zeros(numel(pulse),numel(time));
        if strcmp(class(pulse),'Fitted'), frames = {pulse.width_frames};
        else frames = {pulse.dev_frame}; end
        for j = 1:numel(pulse)
            binary(j,frames{j}) = 1;
        end
    end

    function on = highlight_track(pulse,ID)
        % Parse the ID field and match the correct pulse to highlight
        if strcmp(class(pulse),'Fitted')
            on = ismember([pulse.fitID],ID);
        else
            on = ismember([pulse.trackID],ID);
        end
    end
% ---- End subfunctions ----

end % graph