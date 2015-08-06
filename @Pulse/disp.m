function disp(pulse)
%---- Display overloaded method for Pulse ---

if isempty(pulse), return; end
if numel(pulse) > 1
    display(['Array of Pulse objects of size: [' num2str(size(pulse)) ']']);
    return
end

if pulse.bootstrapped
    
    fprintf('\n')
    display('Simulated pulse from: ')
    display('------ Folder ---------- ')
    display(pulse.input.folder2load)
    fprintf('\n')
    
else
    
    fprintf('\n')
    display('------ Folder ---------- ')
    display(pulse.input.folder2load)
    display(['EmbryoID: ' num2str(pulse.embryoID)])
    display('------ Tracked pulses ---------- ')
    display(['Total tracked pulses: ' num2str(numel(pulse.tracks))])
    display('------ Fitted pulses ----------- ')
    display(['Total fitted pulses: ' num2str(numel(pulse.fits))])
    display('------ Cells ------------------- ')
    display(['Total number of cells: ' num2str(numel(pulse.cells))])
    display(['Total tracked cells: ' ...
        num2str( numel(pulse.cells([pulse.cells.flag_tracked] == 1)) ) ]);
    display(['Total fitted cells: ' ...
        num2str( numel(pulse.cells([pulse.cells.flag_fitted] == 1)) ) ]);
    fprintf('\n')
    
    display('------ Matching ---------------- ')
    if isfield(pulse.categories,'one2one')
        num_one2one = numel( pulse.categories.one2one );
        foo = [pulse.categories.one2one.trackID];
        bar = find_one2one(pulse.map);
        bar = [bar.trackID];
        
        %                 if any( ~ismember(foo,bar)), keyboard; end
    else
        num_one2one = 0;
    end
    display( ['One-to-one matches: ' num2str(num_one2one)] );
    
    if isfield(pulse.categories,'merge')
        num_merge = numel( pulse.categories.merge );
    else
        num_merge = 0;
    end
    display( ['Merged (by fit): ' num2str(num_merge)] )
    
    if isfield(pulse.categories,'split')
        num_split = numel( pulse.categories.split );
    else
        num_split = 0;
    end
    display(['Split (by fit): ' num2str(num_split)])
    
    if isfield(pulse.categories,'miss')
        num_miss = numel( pulse.categories.miss );
    else
        num_miss = 0;
    end
    display(['Missed (by fit): ' num2str(num_miss)])
    
    if isfield(pulse.categories,'add')
        num_add = numel( pulse.categories.add );
    else
        num_add = 0;
    end
    display(['Added (by fit): ' num2str(num_add)])
    fprintf('\n')
    
end
    
end % display
