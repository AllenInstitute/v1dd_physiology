# v1dd_physiology  
This python package allows the interaction with the Allen Institute V1 DeepDive Physiology data set. This data set contains dense calcium imaging (two-photon and three-photon) of excitatory neurons at the center of mouse primary visual cortex while displaying various visual stimuli. In total, data were collected from four mice expressing GCaMP6s as calcium indicator in excitatory neurons. The base data format loosely follows the [Neurodata Without Borders (nwb)](https://www.nwb.org/) standard.

Importantly, the imaged tissue of one mouse were processed by volumetric EM processing and reconstruction to reveal detailed anatomical structures and synaptic-level connectome.  

The description of the data set can be found [here](https://github.com/zhuangjun1981/v1dd_physiology/blob/main/v1dd_physiology/meta/database_description.md).  

The project white paper can be found [here](https://github.com/zhuangjun1981/v1dd_physiology/blob/main/v1dd_physiology/meta/V1DD_WhitePaper_v6.pdf) and additional meta data can be found [here](https://github.com/zhuangjun1981/v1dd_physiology/tree/main/v1dd_physiology/nwb_building/meta_lims).  

The example for using this package can be found in these [ipython notebooks](https://github.com/zhuangjun1981/v1dd_physiology/tree/main/v1dd_physiology/example_notebooks).  

A recorded video introducing this database can be found [here](https://alleninstitute.sharepoint.com/sites/TheBrainObservatoryAnalysisSuperTeam/Shared%20Documents/General/Recordings/Group%20Meeting-20220701_120405-Meeting%20Recording.mp4?web=1) or as a mp4 file: `\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\v1dd_physiology_intro_JunZhuang-20220701.mp4`

&nbsp;  
## Special dependency  
* [NeuroAnalysisTools](https://github.com/zhuangjun1981/NeuroAnalysisTools)  

&nbsp;  
## Install  
#### Create environment  
```
> conda create -n v1dd python=3
> conda activate v1dd
```

#### Install volume_imaging_2P_analysis  
```
> pip install v1dd_physiology
```  
  
or  

```
> git clone https://github.com/AllenInstitute/v1dd_physiology.git
> cd v1dd_physiology
> pip install .
``` 
  
#### Install Jupyter Notebook
```
> conda install jupyter jupyterlab
> python -m ipykernel install --user --name v1dd --display-name "Python3 (v1dd)"
```