DCEMRI.jl
=========

A Fast, Validated Toolkit for Dynamic Contrast Enhanced MRI Analysis

## Why Julia?

From the [Julia website](http://julialang.org/),

> Julia is a high-level, high-performance dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments. It provides a sophisticated compiler, distributed parallel execution, numerical accuracy, and an extensive mathematical function library. The library, largely written in Julia itself, also integrates mature, best-of-breed C and Fortran libraries for linear algebra, random number generation, signal processing, and string processing.

Put simply, it looks like Matlab, which is simple to learn and familiar to most MRI researchers, but it works better and faster and is completely free.  In particular, for the problem of DCE MRI, Julia's simple and flexible parallel computing model allows almost perfect parallelization of the nonlinear least squares fitting problem.  In my informal testing, the intrinsic speed of Julia coupled to my parallel implementation produced a factor of 20-40 speedup over comparable Matlab and Python.

## Installation

### Julia

I've tried to keep the software dependencies to a minimum.  So to run this code you must install only the Julia programming language. Julia can be downloaded from [julialang.org](http://julialang.org/downloads/) as pre-compiled binaries or cloned from the [github repository](https://github.com/JuliaLang/julia).  If you clone the Julia github repo, you should be able to make it with a simple `make` command in the top-level source directory.  The compilation might take a while, but it is completely automatic.  See the [Julia readme](https://github.com/JuliaLang/julia/blob/master/README.md) for compilation instructions. In particular, you must have a [number of development tools](https://github.com/JuliaLang/julia#required-build-tools-and-external-libraries) installed. 

Once you have the base Julia install working, go ahead and start it up.  You might need to set your PATH shell variable to point to your install location, or you may just be able to double-click the icon if you installed one of the pre-compiled binaries.

### DCEMRI.jl

Next, you'll need to install __DCEMRI.jl__ by entering the following at the Julia prompt:
```
julia> Pkg.clone("http://github.com/davidssmith/DCEMRI.jl")
```
This should grab the latest version of __DCEMRI.jl__ and all of the required dependencies.

Now you're ready to run it!

## A Note about Units

All units in the code are SI where possible.  Sometimes, due to numerical accuracy issues, they have been converted internally. But all data should be supplied to the code in SI units.  In particular, time should be in seconds, and relaxation rates in inverse seconds.  Flip angles should be in degrees. The one exception to this in the output is that in the Tofts-Kety models Ktrans is in min^-1.

## Running the Code

### As a shell command

__DCEMRI.jl__ has two basic modes of operation.  The first is command-line invocation, like an operating system command.  To run it as a shell command, first edit the first line of `bin/dcefit` to point to where you installed your Julia binary, as in

```
#!/path/to/julia/binary


```
Next, make sure `bin/dcefit` is executable.  It should already be, but it doesn't hurt to check. Next copy it to a directory that is in your shell's search path.  A good place on UNIX systems, such as Mac OS X, is `/usr/local/bin`.

`dcefit` can parse arguments passed on the command line to configure the model and point to the input data and output file.  To see the available options, run `dcefit -h` at the terminal prompt, you will get

```
usage: dcefit [-O OUTFILE] [-R RELAXIVITY] [-r TR] [-d DCEFLIP]
              [-t T1FLIP [T1FLIP...]] [-m MODELFLAGS] [-p]
              [-w WORKERS] [-v] [-h] [datafile]

Process DCE-MRI data. Optional arguments can be used to override any
values found in input files. For questions, contact David Smith
<david.smith@gmail.com>. For bug reports and feature requests, file an
issue at http://github.com/davidssmith/DCEMRI.jl

positional arguments:
  datafile              path to MAT file containing DCE and T1 data
                        (default: "input.mat")

optional arguments:
  -O, --outfile OUTFILE
                        path to MAT file to contain the ouput
                        (default: "results.mat")
  -R, --relaxivity RELAXIVITY
                        contrast agent relaxivity (1/s) (type:
                        Float64)
  -r, --TR TR           repetition time (ms) (type: Float64)
  -d, --dceflip DCEFLIP
                        flip angle of DCE data (type: Float64)
  -t, --t1flip T1FLIP [T1FLIP...]
                        flip angle(s) of T1 data (type: Float64)
  -m, --modelflags MODELFLAGS
                        logical OR of models to try (1=plasma only,
                        2=Standard, 3=Extended) (type: Int64, default:
                        7)
  -p, --plotting        plot intermediate results
  -w, --workers WORKERS
                        number of parallel workers to use (one per CPU
                        core is best) (type: Int64, default: 4)
  -v, --verbose         show verbose output
  -h, --help            show this help message and exit
```

To process a DCEMRI data set from the command line, the minimum invocation is
`dcefit /path/to/my/dce/data.mat`.

The input data MAT file must contain the following:
- `aif`: Arterial input function (Cp) as a vector, resampled to the DCE time points.
- `dcedata`: DCE data as a 3-D array (1 time by 2 space dimensions).
- `t`: time vector representing the dcedata samples.
- R1 information as either `R1map` and `S0map`, representing pre-calculated R1 relaxation maps, or as `t1data`, indicating that
a multi-flip scan was performed and must be analyzed.  If `t1data` is supplied, the code also needs `t1flip`, a vector of flip angles (in degrees) for the multi-flip data.

The results will be saved in the current directory as `results.mat`.  You can override the output file name and location with the `--outfile` flag.

### As a Julia module




## Validating the Installation

After installing the Julia and the __DCEMRI__ module, you should run the validations, to make sure the calculations work correctly on your machine.  The easiest way to do this is to run `tests/validateall.jl` with a command such as
```
julia -p 4 validateall.jl
```
The `-p 4` flag starts four worker processes.  If that runs, successfully, it will deposit a number of plots into `tests/q4/results` and `tests/q6/results`.  Examine these plots to make sure that the parameters have been recovered accurately.  You can also check the text output of the scripts to see quantitative measures of parameter accuracy.  An example output is shown here:

```
$ julia -p 4 validate6.jl
Running analysis of QIBA v6 (noisefree) data ...
reading input data
found existing R1 map
computing signal enhancement ratios
converting DCE signal to effective R1
converting effective R1 to tracer tissue concentration Ct
fitting DCE data
attempting Standard Tofts-Kety model
fitting 1321 x 8 points on each of 4 workers
processed 30 voxels in 5.6 s (5.3 vox/s)

saving results to results/results.mat
Plotting results ...
Kt RMSE (%): 0.4190757071962333
Kt max error: 2.1663846245047376
Kt CCC: 0.9999304818092664
ve RMSE (%): 0.1261810893311759
ve max error: 0.5704916764566864
ve CCC: 0.9999993740149458

```

To perform the validation on the Quantitative Imaging Biomarkers Alliance phantoms for yourself from the original DICOMS, you will need to download the DICOMS from [Daniel Barboriak's Lab](https://dblab.duhs.duke.edu/modules/QIBAcontent/index.php?id=1).  Then the scripts in the `q4` and `q6` folders will help you translate the DICOM data to MAT files suitable for input into the Julia code.

I have already done this step for you and included the MAT files.  This also avoids you needing to install Python if you don't have it already.  If you want to install Python and run the scripts to convert the DICOM data to MAT files, then I recommend the [Anaconda](http://continuum.io) Python distribution. It has everything you need for scientific programming with Python.

## Running the In Vivo Demo

If you are not already in the __DCEMRI.jl__ source directory, you can navigate there from within
the Julia environment by pressing `;` to obtain the `shell>` prompt and then navigating there
with shell commands.  After each shell command, you will need to press `;` again to drop back into the `shell>` prompt.  The semicolon acts as
a shell command prefix.  The fancy prompt change acts as a reminder that the next text you enter will be interpreted by the shell.

Once you are in the DCEMRI.jl directory, you can run the in vivo data demo with the command
`include("rundemo.jl")`.  After a few seconds to a few minutes, depending on the speed of your machine, you will see the following output text:

```
julia> include("rundemo.jl")
Processing in vivo data ...
reading input data
found multi-flip data
fitting R1 relaxation rate to multi-flip data
fitting 10 x 4582 points on each of 4 workers
processed 18326 voxels in 0.9 s (20983.0 vox/s)

computing signal enhancement ratios
converting DCE signal to effective R1
converting effective R1 to tracer tissue concentration Ct
fitting Standard Tofts-Kety model to tissue concentration Ct
fitting 25 x 1694 points on each of 4 workers
processed 6774 voxels in 0.9 s (7640.2 vox/s)

saving results to results/results.mat
Plotting results ...
elapsed time: 6.123248086 seconds (951780100 bytes allocated, 11.64% gc time)
```

## Concluding Remarks

If you've made it this far, you are ready to run the DCE analysis on your own data.  Congratulations!  If you have problems or find bugs, please file an issue on the [github repository](http://github.com/davidssmith/DCEMRI.jl) or email us.  If you find ways to make it better, please submit your improvements as well. We hope that this can become a community effort that leads to an outstanding, rock solid, trustworthy tool.
