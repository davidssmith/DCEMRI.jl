# This script validate the linearized tofts models

using MAT
using DCEMRI

test_folder = @__DIR__
phantom_ids = (4, 6) # QIBA Phantom #'s
conditions = (:noiseless, :noisy)

# Create noisy data
DCEMRI.makeQibaNoisy(4)
DCEMRI.makeQibaNoisy(6)

for phantom_id in phantom_ids, condition in conditions
  if condition == :noiseless
    mat_filename = "qiba$(phantom_id).mat"
    results_foldername = "results_linear"
  elseif condition == :noisy
    mat_filename = "qiba$(phantom_id)noisy.mat"
    results_foldername = "results_linear_noisy"
  else
    error("Unknown condition: $(condition)")
  end
  qiba_file = joinpath(test_folder, "q$(phantom_id)", mat_filename)
  results_folder = joinpath(test_folder, "q$(phantom_id)", results_foldername)
  results_file = joinpath(results_folder, "results.mat")
  if !isdir(results_folder)
    mkdir(results_folder)
  end
  # The phantom contains 10x10 blocks. Noiseless phantom gobbles them into 1 block.
  voxel_reduction_factor = (condition == :noiseless) ? 10 : 1
  # [5] -> Linearized extended tofts model, [4] -> Linearized tofts model
  models = (phantom_id == 4) ? [5] : [4]

  results = fitdata(datafile=qiba_file, outfile=results_file, models=models)
  DCEMRI.analyze(phantom_id, results, results_folder; dx=voxel_reduction_factor)
end
