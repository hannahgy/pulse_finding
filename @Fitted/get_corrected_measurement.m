function M = get_corrected_measurement(fits,c,meas,input)
%GET_CORRECTED MEASUREMENT Puts a measurement matrix into the
% .corrected_measurement field (placeholder) of an array of Fitted.
%
% USAGE: M = fits.get_corrected_measurement(cells, measMatrix, input)
%
% INPUT: measMatrix - (Nf x p) matrix of measurement to be interpolated
%        input - input array numbered w.r.t. cells.embryoID
% OUTPUT: M - corrected matrix array


fits = fits.align_fits(c,'measurement',meas);
fits = fits.interpolate_traces('measurement',[input.dt]);
M = cat(1,fits.corrected_measurement);

end