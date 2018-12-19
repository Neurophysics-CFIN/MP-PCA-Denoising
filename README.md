# MP-PCA-Denoising

Matlab implementation of Marchenko Pastur denoising (Veraart et al, NeuroImage 142 (2016) 394–406)

Exploiting data redundancy (PCA) and known random matrix properties 
(Marchenko Pastur eigenvalue distribution) to estimate and partially remove noise.

Denoising implementation by Jonas Olesen, Mark Does and Sune Jespersen for diffusion 
MRI data based on the algorithm presented by Veraart et al. (2016) 142, p
394-406 https://doi.org/10.1016/j.neuroimage.2016.08.016.
Modified to remove mean across voxels (compute principal components of
the covariance not correlation matrix).

Free to use, but please cite Veraart et al., NeuroImage (2016) 142, p
394-406 (https://doi.org/10.1016/j.neuroimage.2016.08.016) and Does et al.,
Magn. Reson. Med. (2018) ("Evaluation of Principal Component Analysis Image Denoising on
Multi-Exponential MRI Relaxometry"), reference will be finalized when
available. 
