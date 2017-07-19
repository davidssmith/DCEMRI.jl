# This script validate the linearized tofts models

using MAT
using DCEMRI

# List of tests to run, in pairs
# First run is noiseless, second run is with noise
nList = collect([6,6,4,4])

# Create noisy data
DCEMRI.makeQibaNoisy(4)
DCEMRI.makeQibaNoisy(6)

for i = 1:length(nList)
  n=nList[i]
  # Move to test directory
  cd(Pkg.dir("DCEMRI/test/q$n"))

  # Define the input and output file/directory
  curDx = 1
  if mod(i,2)==1
    datafile="qiba$n.mat"
    outdir = Pkg.dir("DCEMRI/test/q$n") * "/results"
    curDx = 10
  else
    datafile="qiba$(n)noisy.mat"
    outdir = Pkg.dir("DCEMRI/test/q$n") * "/results_noisy"
  end
  isdir(outdir) || mkdir(outdir)

  models=[4]
  if n==4
    models=[5]
  end

  results = fitdata(datafile=datafile, outfile=outdir * "/results.mat", models=models)

  # Make plots and calculate summary statistics
  DCEMRI.analyze(n, results, outdir; dx=curDx)
end
