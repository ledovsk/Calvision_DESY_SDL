# Calvision_DESY_SDL

A simple ROOT macro that allow to inspect individual waveforms from run 200 (DSB crystal at +90 deg) and to perform a  fit of SDL waveform to estimate Cerenkov and Scintillation components

To plot a trigger waveform for event EVT in data file  ```../data/run_200/outfile_LG.root```  and reconstruct trigger timing
```
>. L DSB.C+
> plot_DRS_trigger(EVT,"../data/run_200/outfile_LG.root")
```


To plot DRS waveform in DRS channel CH for event EVT in data file  ```../data/run_200/outfile_LG.root```
```
> plot_DRS_channel(EVT,CH,"../data/run_200/outfile_LG.root")
```

To make Single Delay Line waveform (SDL) with 2ns delay and fit it 
```
reco_WF_SDL(EVT,CH)
```
This macro has location of the data file hardcoded (modify it). It assumes that this file corresponds to DSB crystal.