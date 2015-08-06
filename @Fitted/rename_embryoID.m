function fits = rename_embryoID(fits,embryoID)
% Rename all Fits into a new embryoID
% Please use from PULSE only to ensure CELL objects are
% similarly renamed
old_embryoID = fits(1).embryoID;
[fits.embryoID] = deal(embryoID);

for i = 1:numel(fits)
    
    fID = fits(i).fitID;
    base = 10.^floor(log10(fID) - log10(old_embryoID));
    fID = fID - old_embryoID*base + embryoID*base;
    fits(i).fitID = fID;
    
%     fits(i).stackID = embryoID*1000 + fits(i).cellID;
    
end
end %rename