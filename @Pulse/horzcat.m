function pulse = horzcat(pulse1,pulse2)
%SAFECAT - overloaded concatenation method for putting together two
% Pulse objects. Will look through embryoIDs to ensure there
% are no colliding ID numbers.

pulse = pulse1;
nPulse = numel(pulse); nPulse2 = numel(pulse2);
embryoIDs1 = [pulse1.embryoID]; embryoIDs2 = [pulse2.embryoID];

% find colliding/not colliding embryoIDs
conflicting_embryoIDs = intersect( embryoIDs1, embryoIDs2 );

if isempty( conflicting_embryoIDs)
    for i = 1:nPulse
        pulse( nPulse + i ) = pulse2(i);
    end
else
    
    for i = 1:numel(conflicting_embryoIDs)
        
        thisPulse = pulse2([pulse2.embryoID] == conflicting_embryoIDs(i));
        % new embryoID is +1 to max embryoIDs
        new_embryoID = max([embryoIDs1,embryoIDs2]) + i;
        % reindex fitIDs
        thisPulse = thisPulse.rename_embryoID(new_embryoID);
        
        pulse( nPulse + i ) = thisPulse;
        
    end
    
end

end % cat