## Preparing the Environment

Since we are using an older version of the stack that requires **CentOS**, we need to run a container.  
But before that, we need to gather all required files inside our working directory.

---
💡 Note: Replace any paths (especially the bind mounts) with your own working directory if needed.
---


### 🔧 Create Working Directories

```bash
job_dir="job077321540_mgp8_pp_tttt_wmlep_Q_1000_3000_5f_84TeV"
mkdir "$job_dir"
cd "$job_dir"
mkdir ../mgp8_pp_tttt_wmlep_Q_1000_3000_5f_84TeV
```
### 📦 Copy Required Files

Set the EOS public URL:
```
export EOS_MGM_URL="root://eospublic.cern.ch"
```
### 🧩 MadGraph File
```
cp /eos/experiment/fcc/hh/generation/lhe/mg_pp_tttt_wmlep_Q_1000_3000_5f_84TeV/events_077321540.lhe.gz .
gunzip -c events_077321540.lhe.gz > events.lhe
```

### 🗂 Delphes Cards
```
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/FCChhTrackCov_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/FCChhTrackCov_II.tcl card.tcl
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/muonMomentumResolution_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/trackMomentumResolution_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/electronMomentumResolution_II.tcl .
cp /eos/experiment/fcc/hh/utils/edm4hep_output_config.tcl .
```


### 🧪 PythiaDelphes Files
```
cp /eos/experiment/fcc/hh/utils/config/PythiaDelphes_config_v02.py config.py
cp /eos/experiment/fcc/hh/utils/pythiacards//p8_pp_default.cmd card.cmd
```


### ✏️ Pythia Configuration
```
echo "Beams:LHEF = events.lhe" >> card.cmd
echo "Random:seed = 77321540" >> card.cmd
echo "Main:numberOfEvents = 1000" >> card.cmd
```


### 🚀 Launch the CentOS Container

This command runs the CentOS 7 container with proper volume bindings:
```

apptainer exec \
  --bind /cvmfs/sw.hsf.org:/cvmfs/sw.hsf.org \
  --bind /afs/cern.ch/work/h/hfatehi/centos:/mnt/centos \
  --bind /usr/include:/usr/include \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/library/centos:centos7 \
  /bin/bash
```

### ✅ Inside the Container
Check you're inside CentOS:
```
cat /etc/os-release
```
Verify access to the Key4HEP stack:

```
cat /cvmfs/sw.hsf.org/key4hep/releases/2023-06-05-fcchh/x86_64-centos7-gcc12.2.0-opt/key4hep-stack/*/setup.sh
```
If the above command works, proceed to run your `gen.sh` script.


### ✏️ gen.sh Configuration

```
#!/bin/bash

# Clean environment variables
unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH

# Load Key4HEP stack
source /cvmfs/sw.hsf.org/key4hep/releases/2023-06-05-fcchh/x86_64-centos7-gcc12.2.0-opt/key4hep-stack/*/setup.sh

# Run Delphes with Pythia and EDM4HEP output
DelphesPythia8_EDM4HEP card.tcl edm4hep_output_config.tcl card.cmd events_077321540.root
```



