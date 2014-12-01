from pylab import *
from scipy.io import loadmat, savemat
import time
import sys
import dicom


Rel = 4.5   # assumed Relaxivity of Gd-DTPA at 3 T [s^-1 [mmol Gd-DTPA]^{-1}]
dceflip = 25.0*pi/180.0 # rad
TR = 5e-3   # sec

removektzero = True
noise_sigma = 0.2
random_seed = 1337
seed(random_seed)
deltat = 1
deltax = 1

data_dir = 'DICOM/'
t1_file_ext = 'fa'
dce_file_ext = 'dcemri_testversion4'
outfile_base = 'qiba4'


##
## Extract T1 mapping data
##
nflip = 6
nx = 200
ny = 60
t1flip = array([3.0, 6.0, 9.0, 15.0, 24.0, 35.0])
t1data_dicom = zeros((nflip, nx, ny))
print 'reading DICOMs from', data_dir
for k in range(nflip):
    file_name = '%s/%s%d.dcm' % (data_dir, t1_file_ext, t1flip[k])
    dcm = dicom.read_file(file_name)
    t1data_dicom[k,:,:] = dcm.pixel_array.astype('float')


t1data = t1data_dicom[:,:180,:]
t1data = t1data[:,::deltax,::deltax]
print 'T1 DICOM shape:', t1data.shape

##
## Extract DCE data
##

nt = 661
t = 0.5*arange(nt) + 0.5
dcedata_dicom = zeros((nt, nx, ny))
for k in range(nt):
    file_name = '%s/%s_%03d.dcm' % (data_dir, dce_file_ext, k+1)
    dcm = dicom.read_file(file_name)
    dcedata_dicom[k,:,:] = dcm.pixel_array.astype('float')

dcedata_dicom = dcedata_dicom[:,::deltax,::deltax]
print 'DCE DICOM shape:', dcedata_dicom.shape


data_dce = dcedata_dicom[:,:(180/deltax),:]
nt, nx, ny = data_dce.shape
noise_sigma *= data_dce[0,0,0]


#S0map = ones((nx, ny)) * 1000.0 #
data_aif = mean(mean(dcedata_dicom[:,(180/deltax):,:], axis=2), axis=1)



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
R1p = 1.0 / T1p
Hct = 0.45
S0 = data_aif[:2].mean()
R1_eff_aif = dce_to_r1eff(data_aif, S0, R1p, TR, dceflip)
Cb = r1eff_to_conc(R1_eff_aif.flatten(), R1p, Rel)
Cp = Cb.flatten() / (1.0 - Hct)


mask = zeros((nx,ny)) == 0

data_dce_reduced = data_dce[:,::10,::10].copy()
t1data_reduced = t1data[:,::10,::10].copy()
mask_reduced = mask[::10,::10].copy()

print 'writing MAT files'
mat = {}
mat["R"] = Rel
mat["TR"] = TR
mat["dcedata"] = data_dce_reduced[:,:,1:]
mat["dceflip"] = dceflip*180.0/pi
mat["t1data"] = t1data_reduced[:,:,1:]
mat["t1flip"] = t1flip
mat["t"] = t
mat["aif"] = Cp
mat['mask'] = mask_reduced[:,1:]
mat['extended'] = True
savemat(outfile_base + '.mat', mat)

data_dce = abs(data_dce + noise_sigma*(randn(nt,nx,ny) + 1j*randn(nt,nx,ny)) / sqrt(2.0))
#data_dce = data_dce + noise_sigma*randn(nt,nx,ny) 
mat["dcedata"] = data_dce[:,:,10:]
mat["t1data"] = t1data[:,:,10:]
mat["mask"] = mask[:,10:]
savemat(outfile_base + 'noisy.mat', mat)
