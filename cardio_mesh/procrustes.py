from typing import Dict, List, Optional, Hashable
import numpy as np

from scipy.spatial import procrustes
from scipy.linalg import orthogonal_procrustes
import logging


def mse(s1, s2=None):
    if s2 is None:
        s2 = np.zeros_like(s1)
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


PointCloud = np.ndarray
PointCloudID = Hashable
Transformations = Dict[Literal["traslation", "rotation"], np.ndarray]
ProcrustesResult = Dict[PointCloudID, Transformations]


def generalised_procrustes(
    point_clouds: np.ndarray,
    ids: Optional[List] = None,
    template_mesh: Optional[np.ndarray] = None,
    scaling: bool = False,
    logger: logging.Logger = logging.getLogger()
) -> Dict[Hashable, Dict[str, np.ndarray]]:
    
    """
    Performs Generalized Procrustes Analysis (GPA) to align a set of point clouds 
    to a common reference frame.

    This function iteratively aligns a collection of point clouds to a common template 
    using Procrustes analysis. It applies translation and rotation to minimize the 
    mean squared error (MSE) between the point clouds and the template. Optionally, 
    it can also apply scaling.

    Args:
        point_clouds (np.array): A NumPy array of shape (N, P, D) where:
            - N is the number of point clouds.
            - P is the number of points per cloud.
            - D is the dimensionality of each point (typically 3 for 3D data).
        ids (List, optional): A list of identifiers corresponding to each point cloud. 
            If None, it defaults to `range(len(point_clouds))`.
        template_mesh (np.array, optional): The initial template mesh. If None, 
            the first point cloud is used as the template.
        scaling (bool, optional): Whether to apply scaling during Procrustes 
            alignment. Defaults to False.
        logger (logging.Logger, optional): Logger for logging progress and debug 
            messages. Defaults to `logging.getLogger()`.

    Returns:
        dict: A dictionary mapping each point cloud ID to its transformation 
        parameters. Each entry contains:
            - "traslation" (np.array): The centroid translation applied to the point cloud.
            - "rotation" (np.array): The rotation matrix applied to align the point cloud.
            - (Optional) "scaling" (float): The scaling factor (if scaling is enabled).

    Raises:
        AssertionError: If inputs do not meet expected conditions.

    Notes:
        - If `scaling=True`, Procrustes alignment from `scipy.spatial.procrustes` 
          is used.
        - If `scaling=False`, an orthogonal Procrustes transformation is applied 
          (`scipy.linalg.orthogonal_procrustes`).
        - The algorithm stops iterating when the relative change in disparity is 
          below 1e-2.

    References:
        - https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html
        - https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html
    """

    # Input validation
    assert isinstance(point_clouds, np.ndarray), "point_clouds must be a NumPy array."
    assert len(point_clouds.shape) == 3, "point_clouds must have shape (N, P, D)."
    
    N, P, D = point_clouds.shape
    if ids is None:
        ids = list(range(N))
    
    assert isinstance(ids, list), "ids must be a list."
    assert len(ids) == N, "Length of ids must match the number of point clouds."
    
    if template_mesh is not None:
        assert isinstance(template_mesh, np.ndarray), "template_mesh must be a NumPy array."
        assert template_mesh.shape == (P, D), "template_mesh must have shape (P, D)."

    point_clouds = point_clouds.copy()
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
            if scaling:
                mtx1, mtx2, _disparity = procrustes(template_mesh, point_clouds[i])
                point_clouds[i] = np.array(mtx2)
            else:
                R, s = orthogonal_procrustes(point_clouds[i], template_mesh)
                point_clouds[i] = np.dot(point_clouds[i], R)
                _disparity = mse(point_clouds[i], template_mesh)
                disparity.append(_disparity)

                if it_count == 0:
                    transforms[id]["rotation"] = R
                else:
                    transforms[id]["rotation"] = R.dot(transforms[id]["rotation"])

        template_mesh = point_clouds.mean(axis=0)
        disparity = np.array(disparity).mean(axis=0)
        it_count += 1

    logger.info(
        "Generalized Procrustes analysis with scaling performed after %s iterations"
        % it_count
    )

    return transforms


def transform_mesh(point_cloud: np.ndarray, rotation: Union[None, np.array] = None, traslation: Union[None, np.array] = None):
    
    '''
    params:
    - mesh: Numpy array of size M x 3 representing a point cloud (with M being number of vertices)
    - rotation: rotation Matrix
    - translation: translation vector
    
    returns:
      Numpy array of size M x 3 representing the transformed point cloud       
    '''
       
    point_cloud = point_cloud.copy()
    
    if traslation is not None:
        point_cloud = point_cloud - traslation
        
    if rotation is not None:
        centroid = point_cloud.mean(axis=0)
        point_cloud -= centroid
        point_cloud = point_cloud.dot(rotation)
        point_cloud += centroid
        
    return point_cloud



