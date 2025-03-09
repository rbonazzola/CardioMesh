import os, sys; 

sys.path.append(THIS_DIR := os.path.dirname(os.path.realpath(__file__)))
BASEDIR = os.path.dirname(THIS_DIR)

from CardiacMesh import Cardiac3DMesh
import procrustes as procrustes

import paths as paths
from paths import MESHES_DIR

import Constants

def load_full_heart_mesh(id, timeframe) -> Cardiac3DMesh:
    
    fhm_mesh = Cardiac3DMesh(
        filename = paths.get_3d_pointcloud_file(id, timeframe),
        faces_filename = paths.FACES_FILE,
        subpart_id_filename = paths.SUBPART_ID_FILE)

    return fhm_mesh


def close_chamber(chamber):
    closed_chamber = chamber if chamber.endswith('_closed') else (chamber + "_closed")
    assert closed_chamber in Constants.closed_partitions, f"""
        {closed_chamber} is not cardio_mesh.Constants.closed_partitions. Valid keys are {list(Constants.closed_partitions.keys())}.
    """
    return Constants.closed_partitions[closed_chamber]


def list_mesh_ids():
    return [ file for file in os.listdir(MESHES_DIR) if os.path.isdir(f"{MESHES_DIR}/{file}") ]


def get_aha_thicknesses_for_id(subject_id):
    file = f"{BASEDIR}/data/transforms/thicknesses/{subject_id}_thickness_per_aha.npy"
    assert os.path.exists(file), f"File {file} does not exist."
    return np.load(file)