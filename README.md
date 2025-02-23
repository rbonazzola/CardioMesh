
# CardioMesh
This repository provides a helper class for working with the cardiac meshes produced by the CardioX pipeline.
This code was tested successfully on bi-ventricular meshes produced by the SpASM method as well as the current full-heart model.

# Table of Contents
- [Requirements](#requirements)
- [Tips](#tips)
- [Examples of usage](#Examples-of-usage)
    - [`Cardiac3DMesh`](#cardiac3dmesh)
    - [`Cardiac4DMesh`](#cardiac4dmesh)
    - [`CardiacMeshPopulation`](#cardiacmeshpopulation)
4. [TO DO](#TO-DO)

## Requirements
The code has been tested using the following versions:
- `numpy=1.24.0`,
- `meshio=5.3.5`,
- `vtk=9.3.1`,
- `vedo=2024.5.2`,
- `scipy=1.10.1`,
- `trimesh=4.5.3`
- `pyglet=1.5.18`
- `numpy-stl==2.17.1` (for working with volumetric meshes)


## Installation
While standing in the root directory of this repository, run `pip install -e .`

## Basic description
This repo implements three classes:
- `Cardiac3DMesh`: represents a single cardiac mesh or point cloud.
- `Cardiac4DMesh`: represents a collection of cardiac meshes for one individual, across the cardiac cycle.
- `CardiacMeshPopulation`: represents a population of cardiac meshes (either 3D or 4D), i.e. meshes for different individuals in a population.

## Examples of usage

### `Cardiac3DMesh`

#### Loading a mesh

```
from CardioMesh.CardiacMesh import Cardiac3DMesh

vtk_path = "full_heart_model.vtk"
mesh = Cardiac3DMesh(vtk_path)
mesh.show()
```
![imagen](https://user-images.githubusercontent.com/11581216/124265436-92553100-db2d-11eb-97e0-4227295f1c90.png)

#### Extracting a partition

```
# This follows from the previous code snippet
# LV mesh
lv_mesh = mesh["LV"]

# LV & RV mesh
lvrv_mesh = mesh["LV", "RV"]
```

##### Displaying a mesh (in a Jupyter Notebook)
```
lvrv_mesh.show()
```

![full_heart_model_lvrv](https://user-images.githubusercontent.com/11581216/124301229-6babf000-db57-11eb-8a39-7b3305ae9d89.png)
