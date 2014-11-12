
function dcefit(infile, outfile)
  if nargin < 2, outfile = 'output.mat'; end
  if nargin < 1, infile = 'input.mat'; end
  system(sprintf('dcefit -O %s %s', outfile, infile))
end
