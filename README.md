DCEMRI.jl
=========

A Fast, Validated Open Source Toolkit for Dynamic Contrast Enhanced MRI Analysis

## Why Julia?

From the [Julia website](http://julialang.org/),

> Julia is a high-level, high-performance dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments. It provides a sophisticated compiler, distributed parallel execution, numerical accuracy, and an extensive mathematical function library. The library, largely written in Julia itself, also integrates mature, best-of-breed C and Fortran libraries for linear algebra, random number generation, signal processing, and string processing.

Put simply, it looks like Matlab, which is simple to learn and familiar to most MRI researchers, but it works better and faster and is completely free.  In particular, for the problem of DCE MRI, Julia's simple and flexible parallel computing model allows almost perfect parallelization of the nonlinear least squares fitting problem.  In my informal testing, the intrinsic speed of Julia coupled to my parallel implementation produced a factor of 20-40 speedup over comparable Matlab and Python.

## Installation

Installation is simple.  First you need Julia.
The simplest way to get Julia is to grab the current release version from [julialang.org](http://julialang.org/downloads/).

Next you need __DCEMRI.jl__.  Open Julia. You should see a terminal window with the `julia>` prompt.  This is analogous to the command line in Matlab.
To install __DCEMRI.jl__, run
```
julia> Pkg.clone("http://github.com/davidssmith/DCEMRI.jl")
```
at the `julia>` prompt.  This might take a minute, because the validation data must be downloaded, as well as a few supporting packages. If you get an error than `xcrun` is missing on Mac OS X, follow the instructions to installed the Developer Tools Command Line Tools, and then try again. 

(Optional) Finally, if you want to __DCEMRI.jl__ to create plots for you, you also need Python with [Matplotlib](http://matplotlib.org/) installed.
Most machines probably already have a version of Python with Matplotlib installed.
If you don't have Python with Matplotlib, you can grab the excellent [Anaconda distribution](https://store.continuum.io/cshop/anaconda/), which comes with Matplotlib pre-installed.

Now you're ready to roll!

## A Note about Units

All units in the code are SI where possible.  Sometimes, due to numerical accuracy issues, they have been converted internally. But all data should be supplied to the code in SI units.  In particular, time should be in seconds, and relaxation rates in inverse seconds.  The two exceptions to this rule are that flip angles should be in degrees and Ktrans is output in min^-1.

## Running the Code

### As a Julia module

In the simplest incarnation, if you already have a MAT file containing your data,  you can run the analysis from within Julia using
```
julia> using DCEMRI

julia> results = fitdata(datafile="/path/to/your/datafile.mat")
```
__DCEMRI.jl__ will look for parameters in the input MAT file, and if they are found will use them.  Anything not found in the MAT file will be initialized from the defaults.  These defaults can be viewed with the `defaultparams()` command.  You may also override both the MAT file and the defaults by passing keyword arguments to `fitdata`.

### As a shell command

__DCEMRI.jl__ has two basic modes of operation.  The first is command-line invocation, like an operating system command.  To run it as a shell command, first edit the first line of `bin/dcefit` to point to where you installed your Julia binary, as in

```
#!/path/to/julia/binary


```
Next, make sure `bin/dcefit` is executable.  It should already be, but it doesn't hurt to check. Next copy it to a directory that is in your shell's search path.  A good place on UNIX systems, such as Mac OS X, is `/usr/local/bin`.

`dcefit` can parse arguments passed on the command line to configure the model and point to the input data and output file.  To see the available options, run `dcefit -h` at the terminal prompt, you will get

```
usage: dcefit [-O OUTFILE] [-R RELAXIVITY] [-r TR] [-d DCEFLIP]
              [-c SERCUTOFF] [-t T1FLIP [T1FLIP...]]
              [-m MODELS [MODELS...]] [-p] [-w WORKERS] [-v] [-h]
              [datafile]

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
  -d, --DCEflip DCEFLIP
                        flip angle of DCE data (type: Float64)
  -c, --SERcutoff SERCUTOFF
                        minimum SER to include in processing mask
                        (type: Float64)
  -t, --T1flip T1FLIP [T1FLIP...]
                        list of flip angle(s) of T1 data (type:
                        Float64)
  -m, --models MODELS [MODELS...]
                        list of models: 1=plasma only, 2=Standard,
                        3=Extended (type: Int64)
  -p, --plotting        plot intermediate results
  -w, --workers WORKERS
                        number of parallel workers to use (one per CPU
                        core is good) (type: Int64, default: 4)
  -v, --verbose         show verbose output
  -h, --help            show this help message and exit

```

To process a DCEMRI data set from the command line, the minimum invocation is
`dcefit /path/to/my/dce/data.mat`.

The input data MAT file must contain the following:
- `Cp`: Arterial input function as a vector, resampled to the DCE time points.
- `DCEdata`: DCE data as a 3-D array (1 time by 2 space dimensions).
- `DCEflip` : flip angle in deg of DCE data
- `t`: time vector representing the dcedata samples.
- `TR`: repetition time of DCE scan
- R1 information as either `R10` and `S0`, representing pre-calculated R1 relaxation maps, or as `T1data`, indicating that
a multi-flip scan was performed and must be analyzed.  If `T1data` is supplied, the code also needs `T1flip`, a vector of flip angles (in degrees) for the multi-flip data.

The results will be saved in the current directory as `results.mat`.  You can override the output file name and location with the `--outfile` flag.





## Validating the Installation

After installing the Julia and the __DCEMRI__ module, you should run the validations, to make sure the calculations work correctly on your machine.  The easiest way to do this is to start Julia and then run
```
julia> using DCEMRI

julia> validate()
```
This will run both validations (4 and 6), which could take up to an hour, depending on the number of cores you started Julia with. Examine the results to make sure that the parameters have been recovered accurately.  You can also check the text output of the scripts to see quantitative measures of parameter accuracy.  An example output is shown here:

```
julia> validate(4)
Running analysis of noise-free QIBA v4 data ...
running models
found multi-flip data
fitting R1 relaxation rate to multi-flip data
fitting 6 x 23 points on each of 4 workers
processed 90 voxels in 2.2 s (41.5 vox/s)

computing signal enhancement ratios
converting DCE signal to effective R1
converting effective R1 to tracer tissue concentration Ct
fitting DCE data
attempting Extended Tofts-Kety model
fitting 661 x 23 points on each of 4 workers
processed 90 voxels in 3.8 s (23.5 vox/s)

saving results to /Users/dss/.julia/v0.3/DCEMRI/test/q4/results/results.mat
Plotting results ...
Kt RMSE (%): 6.97465437361441
Kt max error (%): 23.493640353851994
Kt CCC: 0.9998009845162595
ve RMSE (%): 18.02170557638968
ve max error (%): 99.99999999999996
ve CCC: 0.8904290685710147
vp RMSE (%): 23.770196145538407
vp max error (%): 92.10583127104924
vp CCC: 0.9999200988268792
Running analysis of noisy QIBA v4 data ...
running models
found multi-flip data
fitting R1 relaxation rate to multi-flip data
fitting 6 x 2250 points on each of 4 workers
processed 9000 voxels in 0.5 s (19436.3 vox/s)

computing signal enhancement ratios
converting DCE signal to effective R1
converting effective R1 to tracer tissue concentration Ct
fitting DCE data
attempting Extended Tofts-Kety model
fitting 661 x 2250 points on each of 4 workers
processed 9000 voxels in 341.7 s (26.3 vox/s)

saving results to /Users/dss/.julia/v0.3/DCEMRI/test/q4/results_noisy/results.mat
Plotting results ...
Kt RMSE (%): 11.311615941962662
Kt max error (%): 100.0
Kt CCC: 0.9742179876687028
ve RMSE (%): 18.238054961776477
ve max error (%): 100.0
ve CCC: 0.7026132423939505
vp RMSE (%): 12.654024477709797
vp max error (%): 100.0
vp CCC: 0.9717255972607232
Validation complete. Results can be found in /Users/dss/.julia/v0.3/DCEMRI/test/q4.
```

To perform the validation on the Quantitative Imaging Biomarkers Alliance phantoms for yourself from the original DICOMS, you will need to download the DICOMS from [Daniel Barboriak's Lab](https://dblab.duhs.duke.edu/modules/QIBAcontent/index.php?id=1).  Then the scripts in the `q4` and `q6` folders will help you translate the DICOM data to MAT files suitable for input into the Julia code.

I have already done this step for you and included the MAT files.  This also avoids you needing to install Python if you don't have it already.  If you want to install Python and run the scripts to convert the DICOM data to MAT files, then I recommend the [Anaconda](http://continuum.io) Python distribution. It has everything you need for scientific programming with Python.

## Running the In Vivo Demo

You can run the in vivo data demo with the command
`demo()`.  After a few seconds to a few minutes, depending on the speed of your machine, you will see the following output text:

```
julia> demo()
Processing in vivo data ...
running models
found multi-flip data
fitting R1 relaxation rate to multi-flip data
fitting 10 x 4582 points on each of 4 workers
processed 18327 voxels in 1.0 s (19055.6 vox/s)

computing signal enhancement ratios
converting DCE signal to effective R1
converting effective R1 to tracer tissue concentration Ct
fitting DCE data
attempting Standard Tofts-Kety model
fitting 25 x 1694 points on each of 4 workers
processed 6774 voxels in 1.1 s (5928.8 vox/s)

saving results to results/results.mat
Plotting results ...
Demo run complete.
Results can be found in /Users/dss/.julia/v0.3/DCEMRI/demo/results
```

## Concluding Remarks

If you've made it this far, you are ready to run the DCE analysis on your own data.  Congratulations!  If you have problems or find bugs, please file an issue on the [github repository](http://github.com/davidssmith/DCEMRI.jl) or email us.  If you find ways to make it better, please submit your improvements as well. We hope that this can become a community effort that leads to an outstanding, rock solid, trustworthy tool.

To keep your installation of __DCEMRI.jl__ up to date, periodically run `Pkg.update()` at the `julia>` prompt.

## Funding Sources

This project was funded by the National Cancer Institute of the National Institutes of Health, under grants K25CA176219, U01CA142565, R01CA129961, R25CA092043. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
