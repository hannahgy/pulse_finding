function A = get_cell_measurement(pulse,name)
% Returns the specified property of CellObj contained in a pulse (or array
% of Pulses) as a matrix.
%
% USAGE: A = pulses.get_cell_measurement(name)

cells = [pulse.cells];
if isrow(cells(1).(name))
    dim = 1;
else
    dim = 2;
end
A = cat(dim,cells.(name));

end