# Instructions to use for BNL HTCondor farm:

1. Login to eic node
2. Run this command with different parameters (it will use default if not specified)
```bash
TILE_SIZE=12 N_LAYERS=11 SCINTILLATOR_THICKNESS=0.8 ABSORBER_THICKNESS=3 MOMENTUM=1 ./run_submit.sh
```
3. You could modify DEFAULT parameters inside the script.
4. Example in loop:
``` bash
tile_sizes=(10 12 15)
for size in "${tile_sizes[@]}"; do
    export TILE_SIZE=$size
    ./run_submit.sh
done
```
5. This will set everything for you, including `eic-shell`, `epic` repository and `eicrecon`.
And run in batch mode full chain:
```bash
ddsim SIMFILE
eicrecon RECOFILE
root -l -n yourMacro.C (RECOFILE, ANAFILE)
```
6. There is a possibility to wait for jobs using another script `condor_control.sh` which monitors your batch jobs and continues when everything is finished( e.g. you could merge output root analized files which contain your own trees/histograms)

7. The script is a bit long (~500 lines) yet it does fully automated procedure in single file. 
