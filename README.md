# FADE_SAME

**FADE and SAME scores for successful aging in memory**

This code belongs to the paper "A comprehensive score reflecting memory-related fMRI activations and deactivations as a potential biomarker for neurocognitive aging" by Soch, Richter et al. (2021), publicly available from *bioRxiv* and currently under review *Human Brain Mapping*. For instructions how to process these data, see below.

- Preprint: https://www.biorxiv.org/content/10.1101/2021.01.16.426666v1
- Data: TBA
- Code: https://github.com/JoramSoch/FADE_SAME


### Requirements

This code was developed and run using the following software:
- Windows 10 Professional 64-bit
- [MATLAB R2018a/R2020a](https://de.mathworks.com/help/matlab/release-notes.html) (Version 9.4/9.8)
- [MATLAB Statistics and ML Toolbox](https://de.mathworks.com/products/statistics.html) (Version 11.7)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Revision 7771 as of 13/01/2020)
- [spm_helper](https://github.com/JoramSoch/spm_helper) package (Version 1.3)
- [ICC](https://de.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc) function (Version 1.3.1.0)


### Instructions

For re-running analyses reported in the paper, you need to perform the following steps:
1. Create a folder on your computer that hereafter is referred to as the "study directory".
2. Download the [beta images](https://neurovault.org/) (URL TBA) and place them into a sub-folder of the study directory called "beta_images".
3. Download the [analysis scripts](https://github.com/JoramSoch/FADE_SAME/archive/main.zip) and place them into a sub-folder of the study directory called "FADE_SAME".
4. Open MATLAB, set your current directory to this sub-folder, edit the study directory [in line 14](https://github.com/JoramSoch/FADE_SAME/blob/main/analyses_FADE_SAME.m#L14) of `analyses_FADE_SAME.m` and run this script.
5. To reproduce results from the paper, adapt and run the scripts `Figure_*.m` and `Table_*.m`.

* For creating Figure S4, you need to download the MATLAB function `ICC.m` from [MATLAB Central](https://de.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc).