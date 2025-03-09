import meshio
import pandas as pd

import scipy
from scipy import sparse as sp

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

import cardio_mesh
from cardio_mesh.paths import (
    get_3d_pointcloud_file, 
    get_4d_pointcloud
)

LV_OPEN_N_POINTS = 4396
BASE, EPI, ENDO = 0, 1, 2
AHA_FILENAME = f"{cardio_mesh.BASEDIR}/data/LV_{LV_OPEN_N_POINTS}_vertices_with_aha_segments.vtk"
AHA_SEGMENTS_NAMES = list(range(1, 18))


def get_epi_endo_labels():

    import meshio
    
    lv_aha_mesh = meshio.read(AHA_FILENAME)
    lv_aha_labels = lv_aha_mesh.point_data['subpartID'].astype(int)
    
    EPIENDO_FILENAME = f"{cardio_mesh.BASEDIR}/data/LV_4396_vertices_with_epi_endo.vtk"
    epiendo_mesh = meshio.read(EPIENDO_FILENAME)
    epi_endo_labels = epiendo_mesh.point_data['subpartID'].astype(int)
    endo_indices = (epi_endo_labels == ENDO)
    epi_aha_indices = { AHA_SEGMENT: (epi_endo_labels == EPI) & (lv_aha_labels == AHA_SEGMENT) for AHA_SEGMENT in AHA_SEGMENTS_NAMES }

    return epi_aha_indices, endo_indices


def get_lv_indices():
    subpart_df = pd.read_csv(f"{cardio_mesh.BASEDIR}/data/subpartIDs_FHM_10pct.txt", header=None)    
    return (subpart_df == "LV")[0].values


def get_lv_subsetting_mtx():
    subpart_df = pd.read_csv(f"{cardio_mesh.BASEDIR}/data/subpartIDs_FHM_10pct.txt", header=None)    
    lv_indices = (subpart_df == "LV")
    col_ind = lv_indices.index[lv_indices[0]].to_list()
    row_ind = list(range(len(col_ind)))
    
    subsetting_mtx = sp.csc_matrix(
        (np.ones(len(col_ind)), (row_ind, col_ind)), 
        shape=(len(col_ind), subpart_df.shape[0])
    )

    return subsetting_mtx


def load_and_subset_meshes(subject_id, subsetting_mtx=None):
       
    lv_meshes = []
    
    for t in range(1, 50+1):
                
        point_cloud = get_3d_pointcloud(subject_id, t)
        
        if subsetting_mtx is not None:
            point_cloud = subsetting_mtx * point_cloud
        
        lv_meshes.append(point_cloud)
        
    return np.array(lv_meshes)


def compute_thickness_per_aha(point_cloud_array):

    assert isinstance(point_cloud_array, np.ndarray), f"The input point cloud array should be a numpy array. Got {type(point_cloud_array)}."
    assert point_cloud_array.shape[-2] == LV_OPEN_N_POINTS, f"The input point cloud array should have {LV_OPEN_N_POINTS} points (corresponding to a left ventricle without the valve surfaces)."

    if len(point_cloud.shape) == 2:
        thickness_per_segment = []

        for segment in AHA_SEGMENTS_NAMES:
        
            epi_aha_mesh = point_cloud_array[lv_epi_aha_indices[segment]]
            endo_mesh    = point_cloud_array[lv_endo_indices]
            
            epi_aha_mesh_reshaped = epi_aha_mesh.reshape(epi_aha_mesh.shape[0], 1, 3)
            endo_mesh_reshaped    = endo_mesh.reshape(1, endo_mesh.shape[0], 3)
            
            distance_pairs = ((epi_aha_mesh_reshaped - endo_mesh_reshaped)**2).sum(2)
            endo_closest   = distance_pairs.argmin(axis=1)
            
            mean_d = np.sqrt(np.array(
                [ distance_pairs[i, endo_closest[i]] for i in range(distance_pairs.shape[0]) ]
            )).mean()
            
            thickness_per_segment.append(mean_d)
        
        return thickness_per_segment
    
    else: 
        
        if len(point_cloud_array.shape) == 3:
            point_cloud_array = np.expand_dims(point_cloud_array, axis=0)

        n_subjects = point_cloud_array.shape[0]
        n_timeframes = point_cloud_array.shape[1]

        assert len(point_cloud_array.shape) == 4, f"The input point cloud array should have 4 dimensions. Got {point_cloud_array.shape}."    
        thickness_per_segments = []
        for i in range(n_subjects):
            thickness_per_segments_for_i = []
            for t in range(n_timeframes):        
                lv_point_cloud = point_cloud_array[i, t]
                thickness_per_segment = compute_thickness_per_aha(lv_point_cloud)

            thickness_per_segments_for_i.append(thickness_per_segment)
        thickness_per_segments.append(thickness_per_segments_for_i)
            
    return thickness_per_segments
    

def get_lv_meshes(id):
    return load_and_subset_meshes(id, subsetting_mtx=get_lv_subsetting_mtx())


lv_epi_aha_indices, lv_endo_indices = get_epi_endo_labels()