function fits = get_fitID(pulse,fitID)
%GET_CELLID Return the Fitted(s) with the given fitID(s) from a single pulse.
% 
% USAGE: cells = pulse.get_fitID(fitID);

assert(numel(pulse) == 1)

fits = pulse.fits;
fits = fits( ismember([fits.fitID],fitID) );

end