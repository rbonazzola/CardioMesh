import os
import sys
import torch
import pytorch3d

import os
import torch
from pytorch3d.io import load_obj, save_obj
from pytorch3d.structures import Meshes
from pytorch3d.utils import ico_sphere
from pytorch3d.ops import sample_points_from_meshes
from pytorch3d.loss import (
  chamfer_distance, 
  mesh_edge_loss, 
  mesh_laplacian_smoothing, 
  mesh_normal_consistency,
)

from torch import Tensor

from tqdm import tqdm

from CardioMesh.CardiacMesh import Cardiac3DMesh

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['savefig.dpi'] = 80
mpl.rcParams['figure.dpi'] = 80

# Set the device
#if torch.cuda.is_available():
# device = torch.device("cuda:0")
#else:
device = torch.device("cpu")
#    print("WARNING: CPU only, this will be slow!")

# os.system("wget https://dl.fbaipublicfiles.com/pytorch3d/data/dolphin/dolphin.obj")


target = Cardiac3DMesh("mallas/lv_source.vtk")
# verts is a FloatTensor of shape (V, 3) where V is the number of vertices in the mesh
# faces is an object which contains the following LongTensors: verts_idx, normals_idx and textures_idx
# For this tutorial, normals and textures are ignored.

# verts = np.load("mallas/lv_target_verts.npy", allow_pickle=True)
# faces_idx = np.load("mallas/lv_target_faces.npy", allow_pickle=True)

trg_obj = "final_target.obj"
verts, faces, aux = load_obj(trg_obj)
faces_idx = faces.verts_idx

# verts = Tensor(verts)
# faces_idx = Tensor(faces_idx)
# faces_idx = faces_idx.to(device)
# verts = verts.to(device)
# We scale normalize and center the target mesh to fit in a sphere of radius 1 centered at (0,0,0). 
# (scale, center) will be used to bring the predicted mesh to its original center and scale
# Note that normalizing the target mesh, speeds up the optimization but is not necessary!
center = verts.mean(0)
verts = verts - center
scale = max(verts.abs().max(0)[0])
verts = verts / scale
# We construct a Meshes structure for the target mesh
src_mesh = Meshes(verts=[verts], faces=[faces_idx])

############################################################

verts = np.load("mallas/lv_source_verts.npy", allow_pickle=True)
faces_idx = np.load("mallas/lv_source_faces.npy", allow_pickle=True)

# trg_obj = "final_target.obj"
# verts, faces, aux = load_obj(trg_obj)
# faces_idx = faces.verts_idx.to(device)
# verts = verts.to(device)

# source = Cardiac3DMesh("mallas/lv_source.vtk")
verts = Tensor(verts)
faces_idx = Tensor(faces_idx)
faces_idx = faces_idx.to(device)
verts = verts.to(device)

# We scale normalize and center the target mesh to fit in a sphere of radius 1 centered at (0,0,0). 
# (scale, center) will be used to bring the predicted mesh to its original center and scale
# Note that normalizing the target mesh, speeds up the optimization but is not necessary!
center = verts.mean(0)
verts = verts - center
scale = max(verts.abs().max(0)[0])
verts = verts / scale
# We construct a Meshes structure for the target mesh
trg_mesh = Meshes(verts=[verts], faces=[faces_idx])

############################################################

# src_mesh = ico_sphere(4, device)

def plot_pointcloud(mesh, filename="plot.png", title=""):
  # Sample points uniformly from the surface of the mesh.
  points = sample_points_from_meshes(mesh, 5000)
  x, y, z = points.clone().detach().cpu().squeeze().unbind(1)    
  fig = plt.figure(figsize=(5, 5))
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter3D(x, z, -y)
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  ax.set_zlabel('y')
  ax.set_title(title)
  ax.view_init(190, 30)
  plt.savefig(filename)
  plt.close()
  # plt.show()

plot_pointcloud(trg_mesh, "target.png", "Target mesh")
plot_pointcloud(src_mesh, "source.png", "Source mesh")

deform_verts = torch.full(src_mesh.verts_packed().shape, 0.0, device=device, requires_grad=True)

optimizer = torch.optim.SGD([deform_verts], lr=0.1, momentum=0.9)

# Number of optimization steps
Niter = 10000
# Weight for the chamfer loss
w_chamfer = 1.0
# Weight for mesh edge loss
w_edge = 0.1 
# Weight for mesh normal consistency
w_normal = 0.5 
# Weight for mesh laplacian smoothing
w_laplacian = 1 
# Plot period for the losses
plot_period = 100
loop = tqdm(range(Niter))

chamfer_losses = []
laplacian_losses = []
edge_losses = []
normal_losses = []

for i in loop:
    # Initialize optimizer
    optimizer.zero_grad()
                
    # Deform the mesh
    new_src_mesh = src_mesh.offset_verts(deform_verts)
                            
    # We sample 5k points from the surface of each mesh 
    sample_trg = sample_points_from_meshes(trg_mesh, 5000)
    sample_src = sample_points_from_meshes(new_src_mesh, 5000)
    # We compare the two sets of pointclouds by computing (a) the chamfer loss
    loss_chamfer, _ = chamfer_distance(sample_trg, sample_src)

    # and (b) the edge length of the predicted mesh
    loss_edge = mesh_edge_loss(new_src_mesh)

    # mesh normal consistency
    loss_normal = mesh_normal_consistency(new_src_mesh)

    # mesh laplacian smoothing
    loss_laplacian = mesh_laplacian_smoothing(new_src_mesh, method="uniform")

    # Weighted sum of the losses
    loss = loss_chamfer * w_chamfer + loss_edge * w_edge + loss_normal * w_normal + loss_laplacian * w_laplacian
                                                                                                        
    # Print the losses
    loop.set_description('total_loss = %.6f' % loss)
    # Save the losses for plotting
    chamfer_losses.append(float(loss_chamfer.detach().cpu()))
    edge_losses.append(float(loss_edge.detach().cpu()))
    normal_losses.append(float(loss_normal.detach().cpu()))
    laplacian_losses.append(float(loss_laplacian.detach().cpu()))
   
    # Plot mesh
    if i % plot_period == 0:
        plot_pointcloud(new_src_mesh, filename="iter_reversed_2_%d.png" % i, title="iter: %d" % i)
   
    # Optimization step
    loss.backward()
    optimizer.step()


# Fetch the verts and faces of the final predicted mesh
final_verts, final_faces = new_src_mesh.get_mesh_verts_faces(0)

# Scale normalize back to the original target size
final_verts = final_verts * scale + center

# Store the predicted mesh using save_obj
final_obj = 'final_target_2.obj'
save_obj(final_obj, final_verts, final_faces)


