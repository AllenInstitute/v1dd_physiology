# v1dd_physiology  
Code for interacting with the Allen Institute V1 DeepDive dataset.  
  
special dependency:
* [NeuroAnalysisTools](https://github.com/zhuangjun1981/NeuroAnalysisTools)  

# Install
### Create environment
```
> conda create -n v1dd python=3
> conda activate v1dd ("source activate v1dd" for mac or linux)
```

### Install volume_imaging_2P_analysis
```
> pip install v1dd_physiology
```  
  
or  

```
> cd ..
> git clone https://github.com/AllenInstitute/v1dd_physiology.git
> cd v1dd_physiology
> python setup.py install
``` 
  
### Install Jupyter Notebook
```
> conda install jupyter jupyterlab
> python -m ipykernel install --user --name v1dd --display-name "Python3 (v1dd)"
```