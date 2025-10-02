library(oro.dicom)
library(oro.nifti)
library(scales)
library(neurobase)
library(fslr)


mridir = file.path("/home/fsluser/Desktop/002_S_0619/")

t1v4_path = file.path(mridir, "Patient_Visit_4.nii.gz")
t1v1_path = file.path(mridir, "Patient_Visit_1.nii.gz")

mridir = file.path("/home/fsluser/Desktop/002_S_0619/")

t1v4_path = file.path(mridir, "Patient_Visit_4.nii.gz")

nii_v4 = readNIfTI(t1v4_path, reorient = FALSE)

orthographic(nii_v4, xyz= c(128,128,180))
bet_v4 = fslr::fslbet(infile = nii_v4, retimg = TRUE)
orthographic(bet_v4, xyz= c(128,128,180))
cog = cog(bet_v4, ceil = TRUE)
cog = paste("-c", paste(cog, collapse= " "))
bet_v4_2 = fslr::fslbet(infile = nii_v4, retimg = TRUE, opts = cog)

orthographic(bet_v4_2, xyz = c(125,130,200))
orthographic(bet_v4_2, xyz = c(125,130,180))

fast_img = fslr::fsl_biascorrect(mp_v1, retimg = TRUE)

fast_v4 = fslr::fast(file = bet_v4_2, outfile = file.path(paste0(mridir, "/Patient4 Data/02_S_0629_BET_v4.nii.gz")))  #change as per patient

ortho2(bet_v4_2, fast ==1, col.y = alpha("red",0.5), text = "SUBJ_CSF_1_v4", xyz = c(128,128,180))
ortho2(bet_v4_2, fast ==2, col.y = alpha("red",0.5), text = "SUBJ_GM_1_v4", xyz = c(128,128,180))
ortho2(bet_v4_2, fast ==3, col.y = alpha("red",0.5), text = "SUBJ_WM_1_v4", xyz = c(128,128,180))
ortho2(bet_v4_2, fast_v4 ==1, col.y = alpha("red",0.5), text = "SUBJ_CSF_1_v4", xyz = c(128,128,180))
ortho2(bet_v4_2, fast_v4 ==2, col.y = alpha("red",0.5), text = "SUBJ_GM_1_v4", xyz = c(128,128,180))
ortho2(bet_v4_2, fast_v4 ==3, col.y = alpha("red",0.5), text = "SUBJ_WM_1_v4", xyz = c(128,128,180))


pve_CSF_v4 = readNIfTI(paste0(mridir, "02_S_0629_BET_v4_pve_0.nii.gz"))
pve_GM_v4 = readNIfTI(paste0(mridir, "02_S_0629_BET_v4_pve_1.nii.gz"))
pve_WM_v4 = readNIfTI(paste0(mridir, "02_S_0629_BET_v4_pve_2.nii.gz"))



pve_CSF_v4 = readNIfTI(paste0(mridir, "/Patient4 Data/02_S_0629_BET_v4_pve_0.nii.gz"))
pve_GM_v4 = readNIfTI(paste0(mridir, "/Patient4 Data/02_S_0629_BET_v4_pve_1.nii.gz"))
pve_WM_v4 = readNIfTI(paste0(mridir, "/Patient4 Data/02_S_0629_BET_v4_pve_2.nii.gz"))


vdim_csf_v4 = prod(voxdim(pve_CSF_v4))
nvoxels_CSF_v4 = sum(pve_CSF_v4>threshold)


vol_pveCSF_v4 = vdim_csf_v4*nvoxels_CSF_v4/1000
threshold = 0.33

vdim_csf_v4 = prod(voxdim(pve_CSF_v4))

nvoxels_CSF_v4 = sum(pve_CSF_v4>threshold)

vol_pveCSF_v4 = vdim_csf_v4*nvoxels_CSF_v4/1000
print(vol_pveCSF_v4)

vdim_gm_v4 = prod(voxdim(pve_GM_v4))

nvoxels_GM_v4 = sum(pve_GM_v4>threshold)
vol_pveGM_v4 = vdim_gm_v4*nvoxels_GM_v4/1000
threshold = 0.33


vdim_wm_v4 = prod(voxdim(pve_WM_v4))
nvoxels_WM_v4 = sum(pve_WM_v4>threshold)
vol_pveWM_v4 = vdim_WM_v4*nvoxels_WM_v4/1000

vdim_wm_v4 = prod(voxdim(pve_WM_v4))
nvoxels_WM_v4 = sum(pve_WM_v4>threshold)


vol_pveWM_v4 = vdim_wm_v4*nvoxels_WM_v4/1000
print(vol_pveCSF_v4)
print(vol_pveGM_v4)
print(vol_pveWM_v4)
mridir = file.path("/home/fsluser/Desktop/002_S_0619/")

t1v1_path = file.path(mridir, "Patient_Visit_1.nii.gz")


nii_v1 = readNIfTI(t1v1_path, reorient = FALSE)
orthographic(nii_v1, xyz= c(128,128,170))

cog = cog(bet, ceil = TRUE)
cog = paste("-c", paste(cog, collapse= " "))

bet_2 = fslr::fslbet(infile = nii_v1, retimg = TRUE, opts = cog)
orthographic(bet_2, xyz = c(125,130,170))

cog = cog(nii_v1, ceil = TRUE)
cog = paste("-c", paste(cog, collapse= " "))

bet_2 = fslr::fslbet(infile = nii_v1, retimg = TRUE, opts = cog)
orthographic(bet_2, xyz = c(128,128,170))
orthographic(bet_v4_2, xyz = c(125,130,180))

dim(nii_v1)
dim(nii_v4)

t1 = neurobase::readnii(t1v1_path)
t1 [t1<0] = 0
ortho2(robust_window(t1))
ortho2(robust_window(t1), xyz = c(128,128,170))
image(robust_window(t1), useRaster = TRUE)
bc_t1 = bias_correct(file = t1, correction = "N4")
ratio = t1/ bc_t1
ortho2(t1, ratio, xyz = c(128,128,170))

ortho2(robust_window(bc_t1), xyz = c(128,128,170))

sub.bias  <- niftiarr(t1, t1-bc_t1)

q = quantile(sub.bias[sub.bias!=0], probs = seq(0,1,by=0.1))

fcol = div_gradient_pal(low = "red", mid = "yellow", high = "blue")

ortho2(t1, bc_t1, col.y = alpha(fcol(seq(0,1,length = 10)),0.5), ybreaks= q, ycolorbar = TRUE, text = paste0("Original Image Minus N4", "\n Bias Corrected Image"), xyz = c(128,1,170))
















fslview(bet_v4_2)

quants = quantile(bet_2[bet_2 > 0], probs = 0.1)

hist(c(bet_2[bet_2 > 0]), breaks = 200)

quants2 = quantile(bet_v4_2[bet_v4_2 > 0], probs = 0.1)

hist(c(bet_v4_2[bet_v4_2> 0]), breaks = 200)

papaya(list(bet_v4_2))

#devtools::install_github("muschellij2/papayar")

ortho2(bet_v4_2, bet_v4_2>600)

bet = fslr::fslbet(infile = nii_v1, retimg = TRUE)
orthographic(bet, xyz= c(128,128,170))
cog = cog(bet, ceil = TRUE)
cog = paste("-c", paste(cog, collapse= " "))
bet_2 = fslr::fslbet(infile = nii_v1, retimg = TRUE, opts = cog)
orthographic(bet_2, xyz = c(125,130,170))
fast = fslr::fast(file = bet_2, outfile = file.path(paste0(mridir, "/Patient1 Data/02_S_0629_BET_.nii.gz")))
fast = fslr::fast(file = bet_2, outfile = file.path(paste0(mridir, "/Patient1 Data/02_S_0629_BET_.nii.gz")))
ortho2(bet_2, fast ==1, col.y = alpha("red",0.5), text = "SUBJ_CSF_1", xyz = c(128,128,170))
ortho2(bet_2, fast ==2, col.y = alpha("blue",0.5), text = "SUBJ_GM_1", xyz = c(128,128,170))
ortho2(bet_2, fast ==3, col.y = alpha("green",0.5), text = "SUBJ_WM_1", xyz = c(128,128,170))
pve_CSF = readNIfTI(paste0(mridir, "/Patient1 Data/02_S_0629_BET__pve_0.nii.gz"))
#pve_GM = readNIfTI(paste0(mridir, "/Patient1_Data/02_S_0629_BET__pve_1.nii.gz"))
#pve_WM = readNIfTI(paste0(mridir, "/Patient1_Data/02_S_0629_BET__pve_2.nii.gz"))
pve_GM = readNIfTI(paste0(mridir, "/Patient1 Data/02_S_0629_BET__pve_1.nii.gz"))
pve_WM = readNIfTI(paste0(mridir, "/Patient1 Data/02_S_0629_BET__pve_2.nii.gz"))

threshold = 0.33
#calculate the product of voxel dimensions
vdim_csf = prod(voxdim(pve_CSF))
#reads in the pvf file for WM
nvoxels_CSF = sum(pve_CSF>threshold)
#Calculate the volume of CSF in mL^3
vol_pveCSF = vdim_csf*nvoxels_CSF/1000


print(vol_pveCSF)
threshold = 0.33
#calculate the product of voxel dimensions
vdim_gm = prod(voxdim(pve_GM))
#reads in the pvf file for WM
nvoxels_GM = sum(pve_GM>threshold)
#Calculate the volume of CSF in mL^3
vol_pveGM = vdim_gm*nvoxels_GM/1000
print(vol_pveGM)


threshold = 0.33
#calculate the product of voxel dimensions
vdim_WM = prod(voxdim(pve_WM))

#reads in the pvf file for WM
nvoxels_WM = sum(pve_WM>threshold)

#Calculate the volume of CSF in mL^3
vol_pveWM = vdim_WM*nvoxels_WM/1000
print(vol_pveWM)
print(vol_pveCSF)
vdim_WM = prod(voxdim(pve_WM))

#reads in the pvf file for WM
nvoxels_WM = sum(pve_WM>threshold)

#Calculate the volume of CSF in mL^3
vdim_WM = prod(voxdim(pve_WM))

#reads in the pvf file for WM
nvoxels_WM = sum(pve_WM>threshold)

#Calculate the volume of CSF in mL^3
vdim_WM = prod(voxdim(pve_WM))

#reads in the pvf file for WM
nvoxels_WM = sum(pve_WM>threshold)

#Calculate the volume of CSF in mL^3
print(vol_pveCSF)
print(vol_pveGM)
print(vol_pveWM)
