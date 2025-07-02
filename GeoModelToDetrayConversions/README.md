<h1 align="center">
  <br>
  <a href="https://acts.readthedocs.io/en/latest/"><img src="https://avatars.githubusercontent.com/u/48513465?s=200&v=4" alt="ACTS" width="200"></a>
  <br>
  How to make your own ITk geometry and digitization files using a custom ACTS version in Athena
  <br>
</h1>


### Building a custom ACTS version to use in Athena
```bash
$ mkdir athena-custom-acts && cd athena-custom-acts
$ setupATLAS
$ asetup Athena,main,latest
```
> **Note** I used Athena,25.0.34 (2025-06-04T2100 nightly)

```bash
$ lsetup cmake
$ lsetup "views LCG_107a_ATLAS_9 x86_64-el9-gcc13-opt" || true
```

> **Note** LCG version depends on the Athena version! Check in Projects/Athena/build_externals.sh, ```“LCG_${LCG_VERSION_NUMBER}_${LCG_VERSION_POSTFIX}”```
```bash
$ git clone https://github.com/acts-project/acts
```
> **Note** I used commit with hash 0eaaefd2c
```bash
$ cmake -S acts -B acts-build -DCMAKE_INSTALL_PREFIX=$PWD/acts-install -DACTS_BUILD_PLUGIN_JSON:BOOL=ON -DACTS_BUILD_FATRAS:BOOL=ON -DACTS_USE_SYSTEM_GEOMODEL:BOOL=ON -DACTS_BUILD_PLUGIN_TRACCC:BOOL=ON -DCMAKE_CXX_STANDARD=20 -DCMAKE_CUDA_STANDARD=20 -DACTS_BUILD_PLUGIN_GEOMODEL:BOOL=ON
$ cmake --build acts-build
$ cmake --install acts-build
```


### Building local Athena with the custom ACTS version
> **Note** This needs to be a nice, fresh, new terminal
```bash
$ cd Athena-custom-acts
$ setupATLAS
$ asetup --restore
$ export CMAKE_PREFIX_PATH="$PWD/acts-install:$CMAKE_PREFIX_PATH"
$ lsetup git
$ git atlas init-workdir https://gitlab.cern.ch/atlas/athena.git
```
> **Note** Best to use your own fork. You also need a full Athena checkout, or at least all the packages that uses ACTS, please see in the package_filters.txt which ones are necessary to re-build
```bash
$ cd athena/
$ git atlas addpkg EFTrackingDataTransfer
```
> **Note** You can add the piece of code anywhere else you want, I thought EFTracking dir  was convinient. Make sure to change the package_filters if you place it somehwere else. Now place the directory GeoModelConversion from this repository into your sparse checkout of Athena into the athena/Trigger/EFTracking/ directory
```bash
$ cd ..
$ cmake -S athena/Projects/WorkDir -B athena-build -DATLAS_PACKAGE_FILTER_FILE=$PWD/package_filters.txt
```
> **Note** make sure you re-compile all the components that use ACTS in Athena now that you have built a different version to that in the externals, I place here the package_filters I used
```bash
$ source athena-build/x86_64-el9-gcc13-opt/setup.sh
$ export LD_LIBRARY_PATH=$PWD/acts-install/lib64:$LD_LIBRARY_PATH
$ cmake --build athena-build
```

### Dump the files
```bash
$ Reco_tf.py --CA --preInclude "InDetConfig.ConfigurationHelpers.OnlyTrackingPreInclude" --inputRDOFile '/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/PhaseIIUpgrade/RDO/ATLAS-P2-RUN4-03-00-00/mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8481_s4149_r14700/RDO.33629020._000047.pool.root.1' --outputAODFile AOD.test.root --steering doRAWtoALL  --postInclude "GeoModelConversion.GeoModelConversionConfig.GeoModelConversionAlgCfg" --maxEvents 1
```


## Key Dependancies

* [ACTS](https://github.com/acts-project)
  - Documentation [here](https://acts.readthedocs.io/en/latest/)
* [Traccc](https://github.com/acts-project/traccc)
* [Athena](https://gitlab.cern.ch/atlas/athena)
  - Tutorials [here](https://atlassoftwaredocs.web.cern.ch/athena/athena-intro/)

## Support

<a href="https://www.linkedin.com/in/nribaric/" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/purple_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>
