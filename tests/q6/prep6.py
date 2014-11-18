from pylab import *
from scipy.io import loadmat, savemat
import time
import dicom


Rel = 4.5   # assumed Relaxivity of Gd-DTPA at 3 T [s^-1 [mmol Gd-DTPA]^{-1}]
flip_angle = 30 * pi / 180.0 # rad
TR = 5e-3   # sec
nx = 80
ny = 50
nt = 1321
noise_sigma = 0.2
random_seed = 1337
seed(random_seed)
dx = 1
deltat = 1

data_dir = 'DICOM/'
file_ext = 'QIBA_v06_Tofts_beta1'
outfile_base = 'qiba6'

data_dicom = zeros((nt, nx, ny))

t = 0.5*arange(nt)  # ms
print 'reading DICOMs from', data_dir
for k in range(nt):
    file_name = '%s/%s_%04d.dcm' % (data_dir, file_ext, k+1)
    dcm = dicom.read_file(file_name)
    data_dicom[k,:,:] = dcm.pixel_array.astype('float')


data_dce = data_dicom[:,10:70,:]
nt, nx, ny = data_dce.shape
T1map = ones((nx, ny)) # s
R1map = 1 / T1map
S0map = ones((nx, ny)) * 50000.0 #
data_aif = mean(mean(data_dicom[:,70:,:], axis=2), axis=1)

noise_sigma *= data_dce[0,0,0]

# subsample data to speed up the run
data_dce = data_dce[::deltat,:,:]
data_aif = data_aif[::deltat]   # TODO: do this better
t = t[::deltat]
nt = len(t)



# ## 2. Derive the AIF ##
# turn Sb into Cp
print 'converting plasma ROI to AIF'

def dce_to_r1eff(S, S0, R1, TR, flip):
  S = S.T
  S0 = S0.T
  A = S.copy() / S0  # normalize by pre-contrast signal
  E0 = exp(-R1 * TR)
  E = (1.0 - A + A*E0 - E0*cos(flip)) /\
  (1.0 - A*cos(flip) + A*E0*cos(flip) - E0*cos(flip))
  R = (-1.0 / TR) * log(E)
  return R.T

def r1eff_to_conc(R1eff, R1map, relaxivity):
    return (R1eff - R1map) / relaxivity

T1p = 1.440
R1p = 1 / T1p
Hct = 0.45
S0 = data_aif[:4].mean()
R1_eff_aif = dce_to_r1eff(data_aif, S0, R1p, TR, flip_angle)
Cb = r1eff_to_conc(R1_eff_aif.flatten(), R1p, Rel)
Cp = Cb.flatten() / (1.0 - Hct)


## 3. Reduce the problem size averaging 10x10 ROIs to single pixels. ##"
nx /= dx
ny /= dx
data_dce = data_dce[:,::dx,::dx]

R1map_reduced = R1map[::10,::10].copy()
S0map_reduced = S0map[::10,::10].copy()
data_dce_reduced = data_dce[:,::10,::10].copy()

mask = zeros_like(R1map) == 0
mask_reduced = mask[::10,::10].copy()

print 'writing MAT files'
mat = {}
mat["R"] = 4.5
mat["TR"] = 5e-3
mat["dcedata"] = data_dce_reduced
mat["dceflip"] = 30.0
mat["R1map"] = R1map_reduced
mat["S0map"] = S0map_reduced
mat["t"] = t
mat["aif"] = Cp
mat['mask'] = mask_reduced
savemat(outfile_base + '.mat', mat)

data_dce = abs(data_dce + noise_sigma*(randn(nt, nx, ny) + 1j*randn(nt, nx, ny)) / sqrt(2.0))
#data_dce = data_dce + noise_sigma*randn(nt, nx, ny) 
mat["R1map"] = R1map
mat["S0map"] = S0map
mat['dcedata'] = data_dce
mat['mask'] = mask
savemat(outfile_base + 'noisy.mat', mat)

