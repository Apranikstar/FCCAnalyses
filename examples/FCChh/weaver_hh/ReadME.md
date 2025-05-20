### Since we are using an older version of the stack that needs CentOS we need to run a container,
### But before that we have to get all of the files we may need in our container


Crating our working directories
job_dir="job077321540_mgp8_pp_tttt_wmlep_Q_1000_3000_5f_84TeV"
mkdir "$job_dir"
cd "$job_dir"
mkdir ../mgp8_pp_tttt_wmlep_Q_1000_3000_5f_84TeV
--------------------------
now copy the required files:
export EOS_MGM_URL="root://eospublic.cern.ch"
### Madgraph file
cp /eos/experiment/fcc/hh/generation/lhe/mg_pp_tttt_wmlep_Q_1000_3000_5f_84TeV/events_077321540.lhe.gz .
gunzip -c events_077321540.lhe.gz > events.lhe
--------------------------
### Getting Delphes cards usually available at `$DELPHES_DIR/cards/FCC/scenarios/FCChh_II.tcl cards`

cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/FCChhTrackCov_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/FCChhTrackCov_II.tcl card.tcl
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/muonMomentumResolution_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/trackMomentumResolution_II.tcl .
cp /cvmfs/sw.hsf.org/key4hep/releases/2025-01-28/x86_64-almalinux9-gcc14.2.0-opt/delphes/3.5.1pre12-e4qfky/cards/FCC/scenarios/electronMomentumResolution_II.tcl .
cp /eos/experiment/fcc/hh/utils/edm4hep_output_config.tcl .
--------------------------


### Time to get the PythiaDelphes files:

cp /eos/experiment/fcc/hh/utils/config/PythiaDelphes_config_v02.py config.py
cp /eos/experiment/fcc/hh/utils/pythiacards//p8_pp_default.cmd card.cmd

### Pythia Configurations

#### Initial beam file from madgraph
echo "Beams:LHEF = events.lhe" >> card.cmd
#### Setting the seed for Pythia random number generator
echo "Random:seed = 77321540" >> card.cmd
#### Number of events to be generated
echo "Main:numberOfEvents = 1000" >> card.cmd
--------------------------

## Now we are ready to get into our container:


This will automatically run everything you need:

remember that:
    1. the first bind just mounts cvmfs as we use it
    2. second bind gives the container a read-write directory to create files inside, so give replace it with your work directory path
    3. then it loads the centos7 container
    4. it runs that container's bash

So we have:
'''    
apptainer exec \
  --bind /cvmfs/sw.hsf.org:/cvmfs/sw.hsf.org \
  --bind /afs/cern.ch/work/h/hfatehi/centos:/mnt/centos \
  --bind /usr/include:/usr/include \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/library/centos:centos7 \
  /bin/bash
'''

now write this to make sure you have centos up and running:
`cat /etc/os-release`

now we want to check if we have access to our previous stack. to do so we'll just make sure it's setup file returns proper instructions ( a couple of pages of Text):

`cat /cvmfs/sw.hsf.org/key4hep/releases/2023-06-05-fcchh/x86_64-centos7-gcc12.2.0-opt/key4hep-stack/*/setup.sh`

if that works, we'll run the gen.sh file.
--------------------------

The gen.sh file must contain the following code:

#!/bin/bash

### clean the enviroment variables
unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH
### Load 2023-06-05 version of Key4HEP stack
source /cvmfs/sw.hsf.org/key4hep/releases/2023-06-05-fcchh/x86_64-centos7-gcc12.2.0-opt/key4hep-stack/*/setup.sh
#### Using K4SimDelphes to generate EDM4HEP files 
    ##### Usage: DelphesPythia8config_file output_config_file pythia_card output_file

DelphesPythia8_EDM4HEP card.tcl edm4hep_output_config.tcl card.cmd events_077321540.root
--------------------------
