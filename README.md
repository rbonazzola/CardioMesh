# VTK Helpers
This repository provides a helper class for working with the VTK cardiac meshes produced by the CardioX pipeline.
This code was tested successfully on bi-ventricular meshes produced by the SpASM method as well as the current DL-based 4-chamber model.

## Requirements
The code has been tested using the following versions:
- `trimesh==3.9.15`
- `vtk==9.0.1`
- `numpy==1.17.4`
- `pyglet==1.5.18`
- `meshio==2.3.0`

## Tips
You can add this repository as a submodule to your main repository, under the `utils` folder (or similar).

## Examples of usage
### Loading a mesh
```
from VTKHelpers.CardiacMesh import Cardiac3DMesh

vtk_path = "full_heart_model.vtk"
mesh = Cardiac3DMesh(vtk_path)
mesh.show()
```
![imagen](https://user-images.githubusercontent.com/11581216/124265436-92553100-db2d-11eb-97e0-4227295f1c90.png)

### Extracting a partition

```
# This follows from the previous code snippet
# LV mesh
lv_mesh = mesh["LV"]
lv_mesh.show()

# LV & RV mesh
lvrv_mesh = mesh["LV", "RV"]
lvrv_mesh.show()
```

![full_heart_model_lvrv](https://user-images.githubusercontent.com/11581216/124301229-6babf000-db57-11eb-8a39-7b3305ae9d89.png)

## TO-DO
You can contribute to this repository by solving some of the following and submitting a pull request:
- The current version of `VTKObject.extractSubpart` does not support saving a VTK file with subpart labels.
- Create the `VTKObject.show` method to avoid the need of explicit instantiation of a `Trimesh` object.
- Eliminate the dependency on the `meshio` library to save the VTK files.
