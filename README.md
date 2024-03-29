# A practical modification to a resting state fMRI protocol for improved characterization of cerebrovascular function (_https://doi.org/10.1016/j.neuroimage.2021.118306_)
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-5-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

## Breathing task code 

The breathing task code can be found in this repo: https://github.com/RayStick/BreathingTasks_PsychoPy. 
The settings used for the breathing task code are listed below.

### Breath_Hold.py

scan_trigger = 5; doRest = 2; tResting = 480 ; trialnum = 3; tPace = 24 ; tBreathPace = 6; tHold = 15; tExhale = 2; tRecover = 4; BH_instructions = 'BREATH-HOLD task \n \nFollow the breathing instructions \n \nBreathe through your nose'; end_exp_key = 'escape'

### Cued_Deep_Breathing.py

scan_trigger = 5; doRest = 2; tResting = 480; trialnum = 2; tStartRest = 28; tGetReady = 2; tCDB = 8; tCDBPace = 4; tFree = 43; CDB_instructions = 'DEEP BREATHING task \n \nTake deep breaths IN and OUT when cued \n \nBreathe through your nose'; end_exp_key = 'escape'

### Fixation.py

RestDuration = 480; end_exp_key = 'escape'

## ANALYSIS CODE

**PLEASE NOTE**: the analysis code in this repository has been shared alongside the peer-reviewed manuscript in order to be transparent with what analysis was performed to produce the results. The code in this repository is not written as a package that is generalizable to different study set-ups, and unfortunately we do not have the permission to share the data. If you use this code, or parts of it, an acknowledgment would be welcome. You are welcome to use and edit it for your data set-up. If you are looking for citable fMRI pre-processing pipelines check out something like fMRIprep (https://github.com/nipreps/fmriprep). For a package that processes physiological data, to prepare them for analysis with fMRI, see the physiopy community github (https://github.com/physiopy). 

### Step 1 - x.PreProc_MRI

Used the default options (except for input and output directories).

### Step 2 - PreProc_PHYS/preprocess_physdata.m

Ran three times for each data segment.

inputfile_txt = physiological text file from the whole scan session \
time_col = 1 \
vol_col = 2 \
sf = 1000 \
TR = 1.2 \
n_TRs = 534 (BH+REST) or 513 (CDB+REST) or 400 (REST) \
prefix = subjectID_BHREST or subjectID_CDBREST or subjectID_REST \
extra_TRs = 17 \
inputfile_json = 0

### Step 3 - PreProc_PHYS/calc_CO2_regressor.m

Ran three times for each data segment.

filename = output from Step 2 \
fs = 1000 \
vol_col = 1 \
CO2_col = 2 \
prefix = subjectID_BHREST or subjectID_CDBREST or subjectID_REST \
HRFconv= 1 \
output_hires = 1 \
td = 0 \
demean = 0 \
peak_check = 1

### Step 4 - PreProc_PHYS/Shift_hires_CO2regressor

fMRI_file = input is the output from Step 1.
CO2_file = input is the output from Step 3.

Ran in five different ways to make five data segments of length 7 mins 50 seconds:

1. Filename in 'BH+REST'; filename out 'BH+REST'; \
cut_TRs_start=10; \
cut_TRs_end=134;

2. Filename in 'CDB+REST'; filename out 'CDB+REST'; \
cut_TRs_start=10; \
cut_TRs_end=113;

3. Filename in 'BH+REST'; filename out 'rest_bh'; \
cut_TRs_start=144; \
cut_TRs_end=0;

4. Filename in 'CDB+REST'; filename out 'rest_cdb'; \
cut_TRs_start=123; \
cut_TRs_end=0;

5. Filename in 'rest'; filename out 'rest'; \
cut_TRs_start=10; \
cut_TRs_end=0;

Remaining inputs (for the main analysis): \
TR = 1.2 \
extra_TRs_half = 17 \
sf = 1000 \
nf = 40 \
bulk_shift = 0 \
save_bulk_hires = 0 \
fine_shift = 1 \
save_fine_hires = 0 \
shift_unit_s = 0.3 \
shift_max_s = 15 \
prefix = subjectID_BHREST_BH_CO2

Also ran with bulk_shift=1 and fine_shift=0 to explore the bulk shift (Fig4).

### Step 5  - GLM_CO2shifts/x.3dREMLfit_shifts

Ran for each data segment, using the same start and end time defined in Step 4.

Input 1 = Output directory of PreProc_MRI \
Input 2 = Output file from PreProc_MRI \
Input 3 = Output file from PreProc_MRI \
Input 4 = Output file from PreProc_MRI \
Input 5 = Output file from Shift_hires_CO2regressor \
Input 6 = Output file from Shift_hires_CO2regressor (demeaned with 1d_tool.py) \
Input 7 = 100 \
Input 8 = 0.3 \
Input 9 = 1 \
Inputs 10 and 11 = Same start and end times specified for the 5 segments in Step 4 \
Inputs 11 = 1.2

(x.3dDeconvolve_shifts was not run, but is an alternative to x.3dREMLfit_shifts)

### Step 6 - GLM_CO2shifts/x.Threshold_Maps

Ran the code for each each data segment.

Input 1 = Output directory (Concat_3dREMLfit) from Step 5 \
Input 2 = File prefix used for the output files from Step 5 \
Input 3 = Output file from PreProc_MRI \
Input 4 = Output file from PreProc_MRI \
Input 5 = 5.1280e-04 \
Input 6 = 390

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center"><a href="linkedin.com/in/rstickland-phd"><img src="https://avatars.githubusercontent.com/u/50215726?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Rachael Stickland</b></sub></a><br /><a href="https://github.com/BrightLab-ANVIL/Stickland_NeuroImage_2021/commits?author=RayStick" title="Code">💻</a> <a href="#ideas-RayStick" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/BrightLab-ANVIL/Stickland_NeuroImage_2021/commits?author=RayStick" title="Documentation">📖</a></td>
      <td align="center"><a href="https://github.com/kristinazvolanek"><img src="https://avatars.githubusercontent.com/u/54590158?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Kristina Zvolanek</b></sub></a><br /><a href="#ideas-kristinazvolanek" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/BrightLab-ANVIL/Stickland_NeuroImage_2021/pulls?q=is%3Apr+reviewed-by%3Akristinazvolanek" title="Reviewed Pull Requests">👀</a></td>
      <td align="center"><a href="brightlab.northwestern.edu"><img src="https://avatars.githubusercontent.com/u/32640425?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Molly Bright</b></sub></a><br /><a href="#ideas-BrightMG" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center"><a href="https://github.com/smoia"><img src="https://avatars.githubusercontent.com/u/35300580?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Stefano Moia</b></sub></a><br /><a href="#ideas-smoia" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center"><a href="https://www.bcbl.eu/es/conocenos/equipo/cesar-caballero-gaudes"><img src="https://avatars.githubusercontent.com/u/7611340?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Cesar Caballero Gaudes</b></sub></a><br /><a href="#ideas-CesarCaballeroGaudes" title="Ideas, Planning, & Feedback">🤔</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!