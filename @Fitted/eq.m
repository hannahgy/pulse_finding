function equality = eq(fit1,fit2)
% Equality comparator for FITTED
% right now slow, expecting array in first argument and single object
% in the second (fit2). Will return equal if the width_frame of two
% fits have overlap > 5 (frames).
if numel(fit1) > 1 && numel(fit2) > 1
    error('Cannot handle TWO array inputs.');
end

% equality = cellfun( ...
%     @(x,y) numel(find( ismember( x, y) )), ...
%     {fit1.width_frames},{fit2.width_frames} );

equality = false(1,numel(fit1));
for j = 1:numel(fit1)
    if fit1(j).cellID == fit2.cellID
        % can't use bsxfun because of un-uniform output
        if numel(fit1(j).width_frames( ...
                ismember(fit1(j).width_frames, fit2.width_frames))) > 5
            equality(j) = 1;
        end
    end
end

end %eq