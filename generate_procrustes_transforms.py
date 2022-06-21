from scipy.linalg import orthogonal_procrustes
import pickle as pkl
from typing import Dict, List
from IPython import embed
import logging
import torch
import numpy as np
import os

def mse(s1, s2=None):
    if s2 is None:
        s2 = torch.zeros_like(s1)
    return ((s1-s2)**2).sum(-1).mean(-1)

def get_3d_mesh(ids, root_folder):
    
    for id in ids:
        npy_file = f"{root_folder}/{id}/models/fhm_time001.npy"
        pc  = np.load(npy_file)
        yield id, pc
        

def get_4d_mesh(ids, root_folder, timepoints=list(range(1,51))):
    
    for id in ids:
        for t in timepoints:
            npy_file = f"{root_folder}/{id}/models/fhm_time{str(t).zfill(3)}.npy"
            pc  = np.load(npy_file)
        yield id, pc



def generalisedProcrustes(point_clouds: np.array, ids: List, template_mesh=None, scaling=False, logger=logging.getLogger()):


            logger.info("Performing Procrustes analysis with scaling")
            if template_mesh is None:
                template_mesh = point_clouds[0]

            old_disparity, disparity = 0, 1  # random values
            it_count = 0
            
            transforms = {}            


            centroids = point_clouds.mean(axis=1)
            for i, id in enumerate(ids):
                point_clouds[i] -= centroids[i] 
                transforms[id] = {}
                transforms[id]["traslation"] = centroids[i]


            while abs(old_disparity - disparity) / disparity > 1e-2 and disparity:

                old_disparity = disparity
                disparity = []

                for i, id in enumerate(ids):

                    # Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html
                    if scaling:
                        mtx1, mtx2, _disparity = procrustes(template_mesh, point_clouds[i])
                        point_clouds[i] = np.array(mtx2)  # if self.procrustes_scaling else np.array(mtx1)

                    else:
                        # Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html
                        # Note that the arguments are swapped with respect to the previous @procrustes function
                        R, s = orthogonal_procrustes(point_clouds[i], template_mesh)
                        # Rotate
                        point_clouds[i] = np.dot(point_clouds[i], R)  # * s
                        # Mean point-wise MSE
                        _disparity = mse(point_clouds[i], template_mesh) 
                        disparity.append(_disparity)

                        if it_count == 0:
                            transforms[id]["rotation"] = R #, "scaling": s}
                        else:
                            transforms[id]["rotation"] = R.dot(transforms[id]["rotation"]) #, "scaling": transforms[i]["scaling"] * s}

                template_mesh = point_clouds.mean(axis=0)
                disparity = np.array(disparity).mean(axis=0)
                it_count += 1
                
            #self.procrustes_aligned = True
            logger.info(
                "Generalized Procrustes analysis with scaling performed after %s iterations"
                % it_count
            )

            return transforms


root_folder = "/home/home01/scrb/nobackup/meshes/bvalues/Results"
N = 40000
ids = os.listdir(root_folder)[:N]
timepoints=list(range(1,51))
t = 1

point_clouds = []
valid_indices = []
files = [ f"{root_folder}/{id}/models/LV_time{str(t).zfill(3)}.npy" for id in ids ]

for i, id in enumerate(ids):
  try:
      point_clouds.append(np.load(files[i]))
      valid_indices.append(id)
  except:
      pass

point_clouds = np.array(point_clouds)
transforms = generalisedProcrustes(point_clouds, valid_indices)

with open("procrustes_transforms_35k.pkl", "wb") as ff:
    pkl.dump(transforms, ff)

