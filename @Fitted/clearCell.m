function fit_array = clearCell(fit_array)
% Clear all records of fits and tracks from a cell
for i = 1:numel(fit_array)
    
    fitOI = fit_array(i).copy;
    % reset all attributes relating to cell
    fitOI.category = NaN;
    fitOI.bootstrapped = 1;
    fitOI.cellID = NaN;
%     fitOI.stackID = [];
    
    fit_array(i) = fitOI;
    
end
end