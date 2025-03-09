PKG_ROOT_DIR = None
DATA_DIR = None
CACHE_DIR = None
MESHES_DIR = None
FACES_FILE = None
SUBPART_ID_FILE = None

def _initialize_paths():
    import os

    global PKG_ROOT_DIR, DATA_DIR, CACHE_DIR, MESHES_DIR, FACES_FILE, SUBPART_ID_FILE

    PKG_ROOT_DIR = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    DATA_DIR = os.path.join(PKG_ROOT_DIR, "data")
    CACHE_DIR = os.path.join(DATA_DIR, "cached")

    meshes_dir = os.getenv("MESHES_DIR", f"{DATA_DIR}/meshes")
    assert os.path.exists(meshes_dir), f"""
        The folder {meshes_dir}, where meshes should be stored, does not exist. 
        Set the environment variable MESHES_DIR accordingly.    
    """
    MESHES_DIR = meshes_dir

    FACES_FILE = os.path.join(DATA_DIR, "faces_fhm_10pct_decimation.csv")
    SUBPART_ID_FILE = os.path.join(DATA_DIR, "subpartIDs_FHM_10pct.txt")

_initialize_paths()


def get_subsetting_matrix_file(partition): 
    import os
    file = os.path.join(CACHE_DIR, f"subsetting_matrix_{partition}.pkl")
    assert os.path.exists(file), f"{file} does not exist."
    return file


def get_subsetting_matrix(partition): 
    import pickle as pkl
    return pkl.load(open(get_subsetting_matrix_file(partition), "rb"))


def get_mean_shape_file(partition): 
    import os
    file = os.path.join(CACHE_DIR, f"mean_shape_{partition}.npy")
    assert os.path.exists(file), f"{file} does not exist."
    return file


def get_mean_shape(partition):
    import numpy as np
    return np.load(get_mean_shape_file(partition))


def get_procrustes_file(partition):         
    import os
    file = os.path.join(CACHE_DIR, f"procrustes_transforms_{partition}.pkl")
    assert os.path.exists(file), f"{file} does not exist."
    return file


def get_procrustes_transforms(partition):         
    import pickle as pkl
    return pkl.load(open(get_procrustes_file(partition), "rb"))


def get_3d_pointcloud_file(id, timeframe: int):
    import os
    file = f"{MESHES_DIR}/{id}/models/FHM_res_0.1_time{str(timeframe).zfill(3)}.npy"
    assert os.path.exists(file), f"File {file} does not exist"
    return file


def get_3d_pointcloud(ids, timeframe, return_id=False):    
    import numpy as np
    for id in ids:
        npy_file = get_3d_pointcloud_file(id, timeframe)
        point_cloud = np.load(npy_file)
        if return_id:
            yield id, point_cloud
        else:
            yield point_cloud


def get_4d_pointcloud(ids, timepoints=list(range(1, 51)), return_id=False):

    if isinstance(ids, str):
        ids = [ids]

    import numpy as np

    for id in ids:
        dynamic_point_cloud = np.array([np.load(get_3d_pointcloud_file(id, t)) for t in timepoints])
        if return_id:
            yield id, dynamic_point_cloud
        else:
            yield dynamic_point_cloud