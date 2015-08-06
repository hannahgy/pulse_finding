function cellobj = removeFit(cellobj,fit)
%REMOVEFIT - Removes all references of a given FITTED object from a
%CellObj.
%
% USAGE: cellobj.removeFit( fitID );

%@CellObj.removeFit remove a fitID from a cell
fitID = fit.fitID;
idx = find([cellobj.fitID] == fitID); % find removing index
cellobj.fitID(idx) = []; % remove fitID
cellobj.fit_gausses(:,idx) = []; % remove Gaussian colorized peak
cellobj.num_fits = cellobj.num_fits - 1;
% update parameter list (remove)
cellobj.params(:,idx+1) = [];

params = cellobj.params;
t = cellobj.fit_time;
% update residuals
cellobj.residuals = cellobj.raw - ...
    synthesize_gaussians_withbg(params,t);
% update colorized curves
P = plot_peak_color(params(:,2:end),t); % Get colorized fit
cellobj.fit_colorized = P;

end % removeFit