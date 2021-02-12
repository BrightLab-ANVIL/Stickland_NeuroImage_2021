# Stickland_2021

# A practical modification to a resting state fMRI protocol for improved characterization of cerebrovascular function

## Breathing Task code: https://github.com/RayStick/BreathingTasks_PsychoPy

### Breath_Hold.py

scan_trigger = 5  # value the MRI pulse trigger is read in as
doRest = 2  # 0 = no rest; 1 = rest before BH; 2 = rest after BH; 3= rest before AND after BH
tResting = 480  # duration of resting fixation in seconds
trialnum = 3  # number of BH trial repeats
tPace = 24  # duration of paced breathing in seconds
tBreathPace = 6  # duration of each breath in/out in seconds e.g. 6.0 s would be 3s IN and 3s OUT (tPace / tBreathPace pace needs to be integer )
tHold = 15  # duration of BH in seconds
tExhale = 2  # duration for expelling air after BH in seconds
tRecover = 4  # duration of recovery breaths in seconds
BH_instructions = 'BREATH-HOLD task \n \nFollow the breathing instructions \n \nBreathe through your nose'
end_exp_key = 'escape'  # key to press to end the experiment as it is running

### Cued_Deep_Breathing.py

scan_trigger = 5  # value the MRI pulse trigger is read in as
doRest = 2  # 0 = no rest; 1 = rest before CDB task; 2 = rest after CDB task; 3= rest before AND after CDB task
tResting = 480  # duration of resting fixation in seconds
trialnum = 2  # number of CDB repeats
tStartRest = 28  # duration of rest before first CDB section
tGetReady = 2  # duration of get ready warning before CDB section
tCDB = 8  # duration of CDB
tCDBPace = 4  # duration of each breath in/out in seconds e.g. 6.0 s would be 3s IN and 3s OUT (tCDB / tCDBPace needs to be integer )
tFree = 43  # duration of free breathing in between each CDB section
CDB_instructions = 'DEEP BREATHING task \n \nTake deep breaths IN and OUT when cued \n \nBreathe through your nose'
end_exp_key = 'escape'  # key to press to end the experiment as it is running

### Fixation.py

RestDuration = 480  # duration of the resting block in seconds
end_exp_key = 'escape'  # key to press to end the experiment prematurely

## Analysis Code (in this repo)

(to be updated soon)

