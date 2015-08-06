function fits = getFits(cells,fits)

assert(numel(unique([cells.embryoID])) + numel(unique([fits.embryo])) == 2);
assert(unique([cells.embryoID]) == unique([fits.embryo]));

fits = fits( ismember([cells.cellID],[fits.cellID]) );

end