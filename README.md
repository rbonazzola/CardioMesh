# VTK Helpers
This repository provides a helper class for working with the VTK cardiac meshes produced by the CardioX pipeline.
This code was tested successfully on bi-ventricular meshes produced by the SpASM method as well as the current DL-based 4-chamber model.

# Table of Contents
- [Requirements](##Requirements)
- [Tips](#Tips)
- [Examples of usage](#Examples of usage)
    - [`Cardiac3DMesh`](#`Cardiac3DMesh`)
    - [`Cardiac4DMesh`](#`Cardiac4DMesh`)
    - [`CardiacMeshPopulation`](#`CardiacMeshPopulation`)
4. [TO DO](#TO DO)

##Requirements
The code has been tested using the following versions:
- `trimesh==3.9.15`
- `vtk==9.0.1`
- `numpy==1.17.4`
- `pyglet==1.5.18`
- `meshio==2.3.0`

##Tips
You can add this repository as a submodule to your main repository, under the `utils` folder (or similar).

##Examples of usage
---
###`Cardiac3DMesh`
#####Loading a mesh

```
from VTKHelpers.CardiacMesh import Cardiac3DMesh

vtk_path = "full_heart_model.vtk"
mesh = Cardiac3DMesh(vtk_path)
mesh.show()
```
![imagen](https://user-images.githubusercontent.com/11581216/124265436-92553100-db2d-11eb-97e0-4227295f1c90.png)

##### Extracting a partition

```
# This follows from the previous code snippet
# LV mesh
lv_mesh = mesh["LV"]

# LV & RV mesh
lvrv_mesh = mesh["LV", "RV"]
```

#####Displaying a mesh (in a Jupyter Notebook)
```
lvrv_mesh.show()
```

![full_heart_model_lvrv](https://user-images.githubusercontent.com/11581216/124301229-6babf000-db57-11eb-8a39-7b3305ae9d89.png)

---
###`Cardiac4DMesh`
We will assume a folder structure like the following:
```buildoutcfg
├── ED_time_step.csv
├── Ejection_fraction.csv
├── ES_time_step.csv
├── fileCount.txt
├── longitudinalDisplacementEDES.csv
├── massVector.csv
├── output
│   ├── 2D
│   │   ├── LAX_2CH.nii.gz
│   │   ├── LAX_4CH.nii.gz
│   │   └── SAX.nii.gz
│   └── world2gimias
│       ├── output.001.vtk
│       ├── ...
│       └── output.050.vtk
├── quantification2DirectLVVendo.csv
├── quantification2DirectLVVepi.csv
├── quantification2DirectRVV.csv
├── Quantification_temp.txt
├── SQA.csv
├── Stroke_volume.csv
├── vol_per_segm.csv
├── volumeLA.csv
├── volumeLV.csv
├── volumeRA.csv
├── volumeRV.csv
├── wallThickenningEDES.csv
└── wallThickness_per_segm.csv
```

#####Loading a spatio-temporal mesh (a "4D mesh")
```buildoutcfg
vtk_path = "full_heart_model.vtk"
```

---
###`CardiacMeshPopulation`
Here we are assuming a folder structure like the following:
```buildoutcfg

```

```buildoutcfg
vtk_path = "full_heart_model.vtk"
```

## TO-DO
You can contribute to this repository by solving some of the following and submitting a pull request:
- The current version of `VTKObject.extractSubpart` does not support saving a VTK file with non-numeric subpart labels.
- Eliminate the dependency on the `meshio` library to save the VTK files.
