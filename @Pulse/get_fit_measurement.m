function A = get_fit_measurement(pulse,name)
% Returns the specified property of Fitted contained in a pulse (or array
% of Pulses) as a matrix.
%
% USAGE: A = pulses.get_fit_measurement(name)

fits = [pulse.fits];
if isrow(fits(1).(name))
    dim = 1;
else
    dim = 2;
end
A = cat(dim,fits.(name));

end