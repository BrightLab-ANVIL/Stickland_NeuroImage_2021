# Stickland_2021

# A practical modification to a resting state fMRI protocol for improved characterization of cerebrovascular function

## Breathing Task code: https://github.com/RayStick/BreathingTasks_PsychoPy

### Breath_Hold.py

scan_trigger = 5  \
doRest = 2  \
tResting = 480    \
trialnum = 3    \
tPace = 24    \
tBreathPace = 6    \
tHold = 15    \
tExhale = 2    \
tRecover = 4    \
BH_instructions = 'BREATH-HOLD task \n \nFollow the breathing instructions \n \nBreathe through your nose' \
end_exp_key = 'escape'

### Cued_Deep_Breathing.py

scan_trigger = 5  \
doRest = 2  \
tResting = 480   \
trialnum = 2    \
tStartRest = 28    \
tGetReady = 2   \
tCDB = 8    \
tCDBPace = 4    \
tFree = 43    \
CDB_instructions = 'DEEP BREATHING task \n \nTake deep breaths IN and OUT when cued \n \nBreathe through your nose' \
end_exp_key = 'escape'

### Fixation.py

RestDuration = 480  # duration of the resting block in seconds \
end_exp_key = 'escape'  # key to press to end the experiment prematurely

## Analysis Code (in this repo)

(to be updated soon)

