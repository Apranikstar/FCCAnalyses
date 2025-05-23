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
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/FCChh_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/FCChh_II.tcl card.tcl
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


### 🚀 Time to Run Weaver Scripts

Copy the `addons` directory from `FCCAnalyses` to your weaver location.
fix these two lines:
```
from addons.ONNXRuntime.python.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.python.jetClusteringHelper import ExclusiveJetClusteringHelper
```
This output is a bit buggy, remember NOT to put / before afs because the original script does it for you.
```
fccanalysis run stage1.py --output afs/cern.ch/work/h/hfatehi/centos/weaver/test_Hss.root --files-list /afs/cern.ch/work/h/hfatehi/centos/weaver/events_077321540.root --ncpus 16
```
Stage one weaver is now completed!
---------
Stage 2 will work as planned if you change the original script:
---------
change this:
```
from examples.FCCee.weaver.config import variables_pfcand, variables_jet, flavors
```
to this:
```
from config import variables_pfcand, variables_jet, flavors
```

and increase the `maxn = 500` to `maxn = 5000` 

After that just run this script:
```
python stage2.py test_Hss.root out_Hss.root 0 100
```

Congratz, now you have your jet events ready.
---------

Now how to run `weaver-core` using the provided container:

```
#!/bin/bash

CONTAINER_URI='colorsinglet.sif'
singularity shell --nv --bind /afs/cern.ch/work/u/user/weaver-core:/workspace $CONTAINER_URI

export MASTER_PORT=29500
pip install weaver-core
# or just add to PATH if it's missing

torchrun --standalone --nnodes=1 --nproc_per_node=1 \
  -m weaver.train \
  --data-train /workspace/events_077321540.root \
  --data-config /workspace/example.yaml \
  --network-config /workspace/example_ParticleTransformer.py \
  --model-prefix /workspace/trainings/test4GPUs/baseline_test4GPUs \
  --num-workers 1 \
  --gpus 0 \
  --batch-size 2048 \
  --start-lr 1e-3 \
  --num-epochs 20 \
  --optimizer ranger \
  --fetch-step 0.01 \
  --backend nccl


```









