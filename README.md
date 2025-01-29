<h1 align="center">
  <br>
  <a href="https://acts.readthedocs.io/en/latest/"><img src="https://avatars.githubusercontent.com/u/48513465?s=200&v=4" alt="ACTS" width="200"></a>
  <br>
  Traccc in Athena
  <br>
</h1>

<h4 align="center">A minimal example demonstrator for using Traccc with/in Athena.</h4>

<p align="center">
  <a href="#how-to-use">How To Use</a> •
  <a href="#download">Key Dependancies</a> •
  <a href="#license">Support</a>
</p>

## How To Use Traccc in Athena (CPU and GPU versions)

```bash
# Setup environment
$ setupATLAS
$ lsetup git

# Clone this repository
$ git clone https://gitlab.cern.ch/atlas-tdaq-ph2upgrades/atlas-tdaq-eftracking/traccc-integration/gpu
$ cd gpu/G-200

# Create a build directory and go into the directory
$ mkdir build
$ cd build

# Setup Athena release (> 24.0)
# I used Athena,25.0.24 I recommend using this instead of main.
$ asetup Athena,main,latest

# Install the dependancies and compile the code
$ cmake ../traccc-athena
$ make -j4
# Sourcing breaks some paths, so the temporary fix is:
$ env > envlog.log
$ source */setup.sh
$ export `grep CMAKE_PREFIX_PATH envlog.log`

# Run the reco chain
$ cd ../traccc-athena/run/
$ Reco_tf.py --CA --preInclude "InDetConfig.ConfigurationHelpers.OnlyTrackingPreInclude" --inputRDOFile '/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/PhaseIIUpgrade/RDO/ATLAS-P2-RUN4-01-01-00/mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8481_s4038_r14365/RDO.32520484._000001.pool.root.1' --outputAODFile AOD.test.root --steering doRAWtoALL  --postInclude "EFTracking.TrackingAlgConfig.TrackingAlgCfg" --maxEvents 1

# Make pretty plots
$ runIDTPM.py --inputFileNames AOD.test.root --outputFile test_IDTPM_output --trkAnaCfgFile ./EFTracking_TrkAnaConfig.json

```

> **Some useful information** <br>
> Traccc version is specified in the source file. You should update the tag in case you are interested in a different version. <br> <br>
> You can also use a local direcotry, but it requires a lot of massaging, so I do not recommend. The instructions are below if you are nonetheless keen on going on an adventure.<br> <br>
> The detector geometry, detector digitization, surface grids and Athena-to-Detray map files and paths are currently hardcoded, so you have to change to whichever directory you want to use. <br> <br>
> Please note that you will need to be a member of e-group atlas-tdaq-phase2-EFTracking-developers to access /eos/project/a/atlas-eftracking/GPU/ITk_data/ <br> <br>
> You can also recreate the files by using the code [here](https://gitlab.cern.ch/atlas-tdaq-ph2upgrades/atlas-tdaq-eftracking/traccc-integration/Common/-/tree/main/GeoModelToDetrayConversions?ref_type=heads)<br> <br>
>
>
> This repo for now only includes clustering for ITk pixels, note: <br>
> * the measurements are shifted wrt Athena measuremens, because of: <br>
>   * Lorentz shift not being applied
>   * omega weighting of the measurements not being applied, [see here](https://gitlab.cern.ch/nribaric/athena/-/blob/master/InnerDetector/InDetRecTools/SiClusterizationTool/src/ClusterMakerTool.cxx#L289)

## If you want to use a local build of Traccc

```bash
$ setupATLAS
$ asetup none,gcc13,cmakesetup --cmakeversion=3.29.5
```
> **Note:** <br> 
> Please make sure the LCG versions and nightlies are correct for the Athena release you want to use
```bash
$ lsetup git
$ git clone https://github.com/acts-project/traccc.git
$ mkdir build-traccc && cd build-traccc
$ cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_PREFIX_PATH="/cvmfs/atlas-nightlies.cern.ch/repo/sw/main_Athena_x86_64-el9-gcc13-opt/2024-11-12T2101/AthenaExternals/25.0.22/InstallArea/x86_64-el9-gcc13-opt;/cvmfs/sft.cern.ch/lcg/views/LCG_106a_ATLAS_3/x86_64-el9-gcc13-opt"    -DTRACCC_BUILD_CUDA=TRUE -DTRACCC_USE_SYSTEM_TBB=TRUE -DTRACCC_USE_SYSTEM_EIGEN3=TRUE -DTRACCC_USE_SYSTEM_ACTS=TRUE -DALGEBRA_PLUGINS_USE_SYSTEM_VC=TRUE -DDETRAY_USE_SYSTEM_NLOHMANN=TRUE -DVECMEM_BUILD_CUDA_LIBRARY=TRUE -DTRACCC_USE_ROOT=FALSE -DTRACCC_BUILD_TESTING=FALSE -DCMAKE_INSTALL_PREFIX=../install-traccc ../traccc/
```
> **Note:** <br>
> I used a specific commit from Traccc: ```034eef4``` <br>
> For that you have to use a specific detray commit that plays nicely with traccc: <br>
```-DTRACCC_DETRAY_SOURCE="GIT_REPOSITORY;https://github.com/acts-project/detray.git;GIT_TAG;028de5d3910cbe5037a0a2dd2f60de140e1b09ff"``` <br>

```bash
$ make
$ make install
```

> **Note** <br>
> Now in a brand new, fresh terminal, setup Athena

```bash
$ mkdir athena-build && cd athena-build
$ setupATLAS
$ setup Athena,main,latest
```
> **Note** <br> we had to specify ```Athena,main,latest,gcc13``` because by default the latest nightly comes with gcc14

```bash
$ export CMAKE_PREFIX_PATH=$PWD/../install-traccc:$CMAKE_PREFIX_PATH
$ cmake -DATLAS_USE_SYSTEM_TRACCC=TRUE ../traccc-athena/
$ make
$ env > envlog.log
$ source */setup.sh
$ export `grep CMAKE_PREFIX_PATH envlog.log`
```


## Key Dependancies

* [ACTS](https://github.com/acts-project)
  - Documentation [here](https://acts.readthedocs.io/en/latest/)
* [Traccc](https://github.com/acts-project/traccc)
* [Athena](https://gitlab.cern.ch/atlas/athena)
  - Tutorials [here](https://atlassoftwaredocs.web.cern.ch/athena/athena-intro/)

## Support

<a href="https://www.linkedin.com/in/nribaric/" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/purple_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>

