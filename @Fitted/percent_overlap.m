function [perc,varargout] = percent_overlap(fits,cells)
%PERCENT_OVERLAP Counts the percentage of overlapping frames from the
% consecutive pulse-frames used to analyze each pulse.
%
% USAGE: perc = fits.percent_overlap(cells);
% 		 [perc,counts] = fits.percent_overlap(cells);
%
% xies@mit Oct 2013

count = zeros(size(cat(2,cells.area)));

for i = 1:numel(fits)
    
    this_fit = fits(i);
    I = find([cells.stackID] == this_fit.stackID);
    count( this_fit.margin_frames(3:end-2), I) = ...
        count( this_fit.margin_frames(3:end-2), I) + 1;
    
end

perc = numel(count(count > 1)) / numel(count(count > 0));
if nargout > 1, varargout{1} = count; end

end % percent_overlap