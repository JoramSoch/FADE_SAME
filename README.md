# FADE_SAME

**FADE and SAME scores for successful aging in memory**

This code belongs to the paper "A comprehensive score reflecting memory-related fMRI activations and deactivations as a potential biomarker for neurocognitive aging" by Soch, Richter et al. (2021), published open access in *Human Brain Mapping*. For instructions how to process these data, see below.

- Paper: https://onlinelibrary.wiley.com/doi/10.1002/hbm.25559
- Data (1): https://neurovault.org/collections/EPLZNQAD/
- Data (2): https://github.com/JoramSoch/FADE_SAME_SPM_mat
- Code: https://github.com/JoramSoch/FADE_SAME


### Requirements

This code was developed and run using the following software:
- Windows 10 Professional 64-bit
- [MATLAB R2018a/R2020a](https://de.mathworks.com/help/matlab/release-notes.html) (Version 9.4/9.8)
- [MATLAB Statistics and ML Toolbox](https://de.mathworks.com/products/statistics.html) (Version 11.7)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Revision 7771 as of 13/01/2020)
- [spm_helper](https://github.com/JoramSoch/spm_helper) package (as on GitHub)
- [ICC](https://de.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc) function (Version 1.3.1.0)


### Instructions

For re-running analyses reported in the paper, you need to perform the following steps:
1. Create a folder on your computer that hereafter is referred to as the "study directory".
2. Download the [beta images](https://neurovault.org/collections/EPLZNQAD/) and place them into a sub-folder of the study directory called "beta_images".
3. Download the [SPM.mat files](https://github.com/JoramSoch/FADE_SAME_SPM_mat/archive/main.zip) and place them into the same sub-folder of the study directory.
4. Download the [analysis scripts](https://github.com/JoramSoch/FADE_SAME/archive/main.zip) and place them into a sub-folder of the study directory called "FADE_SAME".
5. Open MATLAB, set your current directory to this sub-folder, edit the study directory [in line 14](https://github.com/JoramSoch/FADE_SAME/blob/main/analyses_FADE_SAME.m#L14) of `analyses_FADE_SAME.m` and run this script.
6. To reproduce results from the paper, adapt and run the scripts `Figure_*.m` and `Table_*.m`.

* When creating Figure S3, there will be an error, as ApoE genotype is not included in the public data share. These data will be made available upon request, please write an e-mail to the [corresponding authors of the paper](https://onlinelibrary.wiley.com/doi/10.1002/hbm.25559).
* For creating Figure S4, you need to download the MATLAB function `ICC.m` from [MATLAB Central](https://de.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc).