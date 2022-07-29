import vtk
import numpy as np
import os
import meshio  # tested with 2.3.0

import sys
from Constants import *
from tqdm import tqdm

import logging
import random

from IPython import embed  # For debugging

import pickle as pkl
from stl import mesh as stlmesh

"""
This module is aimed to simplify the implementation of common tasks on VTK triangular meshes,
that self.points overly convoluted if the usual VTK Python wrapper for C++ is used,
and render the code difficult to follow.
"""


def set_logger(logger):
    return logger if logger is not None else logging.getLogger()


class Cardiac3DMesh:

    """
    This class represents a single cardiac mesh or point cloud.
    """

    def __init__(
        self,
        filename=None,
        subpartIDs=None,
        load_connectivity_flag=True,
        dataset_version=None,
        logger=None,
    ):

        """
        filename: path to a file for a VTK Polydata object.
        subpartIDs: either a label or a list of labels to subset the mesh for (e.g. "LV", "RV", etc.)
        load_connectivity_flag: boolean indicating whether or not to load the mesh connectivity information
        dataset_version: if None, it is inferred automatically based on the subpart IDs. Currently supported values are "LEGACY_2CHAMBER_SPASM" and "FULL_HEART_MODEL_MMF"
        """

        self.points = None
        self._logger = set_logger(logger)

        if filename is not None:
            self._filename = filename

            # TODO: inform when a symbolic link is broken
            if not os.path.exists(self._filename):
                raise FileExistsError("File {} does not exist.".format(self._filename))

            # check if filename extension is .vtk or pickle 
            if self._filename.endswith(".vtk"):

                self._reader = vtk.vtkPolyDataReader()
                self._reader.SetFileName(self._filename)
                self._reader.Update()

                self._load_point_cloud()
                if load_connectivity_flag:
                    self._load_connectivity()
                self._load_partition_ids()
                self._infer_dataset_version()

            elif self._filename.endswith(".pkl"):

                with open(self._filename, "rb") as f:
                    dict = pkl.load(f)
                
                # We can assume vertices and faces are numpy arrays
                self.points = dict['points']
                self.triangles = dict['triangles']
                try:
                    self.subpartID = dict['subpartID']
                except:
                    subpartIDs = None
                self._infer_dataset_version()

        if subpartIDs is not None:
            newMesh = self._extract_subpart(subpartIDs)
            self.__dict__.update(newMesh.__dict__)


    def _load_point_cloud(self):

        """
        :return: numpy array where each element is a triple of (x, y, z)
        coordinates and a set of indices representing the links to that point
        """

        output = self._reader.GetOutput()
        n_points = output.GetNumberOfPoints()
        self.points = np.array([output.GetPoint(i) for i in range(n_points)])

    def _load_connectivity(self, triangles=None, edges=None):

        """
        triangles: if not None, must be a list of 3-tuples containing valid point indices
        edges: if not None, must be a list of 2-tuples containing valid point indices
        if both are None, the connectivity will be read from the VTK file
        """

        if triangles is None and edges is None:
            # TODO: raise an error if the VTK file has not been provided
            output = self._reader.GetOutput()
            self.n_cells = output.GetNumberOfCells()
            self.triangles = [
                [int(output.GetCell(l).GetPointId(k)) for k in (0, 1, 2)]
                for l in range(self.n_cells)
            ]
        elif triangles is not None:
            # TODO: perform better sanity check on the `triangles` argument
            if not isinstance(triangles, list):
                raise TypeError
            self.triangles = triangles
        else:
            self.triangles = self._triangles_from_edges(edges)

        self.triangles = np.array(self.triangles)

    def _infer_dataset_version(self):

        """
        Infer the dataset version. Currently, "LEGACY_2CHAMBER_SPASM" and "FULL_HEART_MODEL_MMF" are supported.
        """

        # TODO: decide if there is a better to infer the partition
        if self.distinct_subparts == {1, 2, 4}:
            self._dataset_version = LEGACY_2CHAMBER_SPASM
            self._subpart_id_mapping = LEGACY_2CH_SUBPART_IDS
        else:
            self._dataset_version = FULL_HEART_MODEL_MMF
            self._subpart_id_mapping = FHM_SUBPART_IDS

    @property
    def edges(self):
        try:
            return self._edges
        except AttributeError:
            self._edges = self._edges_from_triangles(self.triangles)
            return self._edges

    def _edges_from_triangles(self, triangles):
        edges = []
        for x, y, z in triangles:
            edges.extend([(x, y), (y, z), (x, z), (y, x), (z, y), (z, x)])
        return list(set(edges))

    def _triangles_from_edges(self):
        raise NotImplementedError
        # TODO: implement

    # @property
    # def triangles(self):
    #     try:
    #         return self._triangles
    #     except AttributeError:
    #         output = self._reader.GetOutput()
    #         self._triangles = [[int(output.GetCell(j).GetPointId(i)) for i in (0, 1, 2)] for j in range(self.n_cells)]
    #         return np.array(self._triangles)

    @property
    def neighbors_dict(self):
        try:
            return self._neighbors_dict
        except AttributeError:
            self._neighbors_dict = {}
            for edge in self.edges:
                self._neighbors_dict.get(edge[0], []).append(edge[1])
            return self._neighbors_dict

    @property
    def distinct_subparts(self):
        return set(self.subpartID)

    @property
    def v(self):
        """
        This works as an alias, for compatibility with third-party code.
        Probably it is better to use a wrapper (child) class for this.
        """
        return self.points

    @property
    def f(self):
        """
        This works as an alias, for compatibility with third-party code.
        Probably it is better to use a wrapper (child) class for this.
        """
        return self.triangles

    def _load_partition_ids(self):

        """
        Generate a list of the subpart IDs for each of the vertices (i.e. which partition of the mesh they belong to)
        """

        output = self._reader.GetOutput()

        pp = output.GetPointData().GetAbstractArray(0)
        self.subpartID = [pp.GetValue(i) for i in range(self.n_points)]
        # self.subpartID = [int(pp.GetComponent(i, 0)) for i in range(self.n_points)]

    @property
    def n_points(self):
        """
        Number of vertices in the mesh
        """
        return len(self.points)

    @property
    def n_edges(self):
        """
        Number of edges in the mesh (without double-counting)
        """
        return len(self.edges)

    def __len__(self):
        return self.n_points

    def __repr__(self):
        return "Point cloud\n\n {} \n\n.".format(
            self.points.__str__()
        )  # with connectivity\n\n{}".format(self.points.__str__(), self.triangles.__str__())

    def show(self, engine="trimesh"):

        if engine == "trimesh":
            try:
                from trimesh import Trimesh

                # return Trimesh(self.v, self.f).show()
                Trimesh(self.v, self.f).show()
            except ImportError:
                logger.error(
                    "Trimesh package must be installed in order to use this functionality."
                )

    def _extract_subpart(self, ids):
        """
        ids: a label or a list of labels for the subpart/s to be extracted
        :return: a Cardiac3DMesh object representing the subpart to be extracted.
        """

        ids = [ids] if not isinstance(ids, list) else ids
        
        subvtk = Cardiac3DMesh()
        subvtk.points = np.array(
            [self.points[i] for i in range(self.n_points) if self.subpartID[i] in ids]
        )
        subvtk.subpartID = np.array(
            [
                self.subpartID[i]
                for i in range(self.n_points)
                if self.subpartID[i] in ids
            ]
        )

        if self.triangles is not None:
            point_ids = [i for i, id in enumerate(self.subpartID) if id in ids]
            point_ids_set = set(point_ids)
            triangles = [
                tuple(triangle)
                for triangle in self.triangles
                if all([pp in point_ids_set for pp in triangle])
            ]

            id_mapping = {x: i for i, x in enumerate(point_ids)}
            subvtk.triangles = np.array(
                [tuple([id_mapping[x] for x in triangle]) for triangle in triangles]
            )

        return subvtk

    def _map_subpart_ids(self, ids):
        """ """

        ids = list(ids) if isinstance(ids, tuple) else ids
        ids = [ids] if not isinstance(ids, list) else ids

        # 1 -> [1]
        # 1,2 -> [1,2]
        # "LV" -> ["LV"]
        # "LV", "RV"  -> ["LV", "RV"]

        possible_values = list(self.distinct_subparts)
        possible_values += [
            x
            for x in self._subpart_id_mapping
            if all([y in self.distinct_subparts for y in self._subpart_id_mapping[x]])
        ]

        kk = []

        for id in ids:
            for x in self._subpart_id_mapping.get(id, [id]):
                if x in possible_values:
                    kk.append(x)
                else:
                    raise ValueError(
                        "{} is not a valid partition (use {} or combinations thereof)".format(
                            x, ", ".join(sorted([str(x) for x in possible_values]))
                        )
                    )

        return kk

    def __getitem__(self, ids):
        return self._extract_subpart(self._map_subpart_ids(ids))

    @property
    def shape(self):
        return self.points.shape

    @property
    def adj_matrix(self):

        """
        Returns a sparse matrix (of size #verts x #verts) where each nonzero
        element indicates a neighborhood relation. For example, if there is a
        nonzero element in position (15,12), that means vertex 15 is connected
        by an edge to vertex 12.
        """

        try:
            return self._adj_matrix
        except AttributeError:
            from scipy import sparse as sp

            self._adj_matrix = sp.csc_matrix(
                (
                    np.ones(len(self.edges)),
                    ([x[0] for x in self.edges], [x[1] for x in self.edges]),
                )
            )
            return self._adj_matrix

    def adj_matrix_to_edges(self, adj_matrix):
        from scipy import sparse as sp

        non_zero_indices = sp.find(adj_matrix)
        return zip(non_zero_indices[0], non_zero_indices[1])

    # mesh to vtk
    def save_to_vtk(self, filename):
        meshio.write_points_cells(
            filename,
            self.points,
            cells={"triangle": np.array(self.triangles)},
            point_data={"subpartID": self.subpartID},
        )
        
    # mesh to pickle
    def save_to_pkl(self, filename):
        dict = {"points" : self.points, "triangles" : self.triangles, "subpartID" : self.subpartID}

        with open(filename, "wb") as f:
            pkl.dump(dict, f)

    # mesh to stl
    def save_to_stl(self, filename):
        num_triangles = self.triangles.shape[0]
        data = np.zeros(num_triangles, dtype=stlmesh.Mesh.dtype)

        for i in range(num_triangles):
            #I did not know how to use numpy-arrays in this case. This was the major roadblock
            # assign vertex co-ordinates to variables to write into mesh
            v1x, v1y, v1z = self.points[self.triangles[i,0],0], self.points[self.triangles[i,0],1], self.points[self.triangles[i,0],2]
            v2x, v2y, v2z = self.points[self.triangles[i,1],0], self.points[self.triangles[i,1],1], self.points[self.triangles[i,1],2]
            v3x, v3y, v3z = self.points[self.triangles[i,2],0], self.points[self.triangles[i,2],1], self.points[self.triangles[i,2],2]
            
            data["vectors"][i] = np.array([[v1x, v1y, v1z],[v2x, v2y, v2z],[v3x, v3y, v3z]])

        m = stlmesh.Mesh(data)
        m.save(filename)

class Cardiac4DMesh:

    """
    Class representing a collection of cardiac meshes for one individual, across the cardiac cycle.
    Public attributes:
      meshes:
      triangles:
      time_frames: list of integers
      subjectID
      LVEDV, LVESV, LVEF, LVSV, LVM
    """

    def __init__(self, root_folder, time_frames=None, logger=None):

        """
        root_folder: path to the PyCardioX output for the given individual
        time_frames: list|tuple containing 1-50 and/or "ED"/"ES"
        """

        self._root_folder = root_folder
        self._logger = set_logger(logger)

        if time_frames is None:
            self.time_frames = [i + 1 for i in range(50)]
        else:
            self.time_frames = time_frames
        self._time_frames_to_int()

        self._vtk_paths = self._get_vtk_paths()
        self._load_meshes()

    # def _set_logger(self, logger):
    #    return logger if logger is not None else logging.getLogger()

    @property
    def subjectID(self):
        """
        This method assumes that the subject ID is the name of the rightmost folder.
        """
        self._subjectID = os.path.basename(self._root_folder.strip("/"))
        return self._subjectID

    def _load_meshes(self, load_connectivity_flag=True):
        self.meshes = []
        for i, vtk_path in enumerate(self._vtk_paths):
            if i == 0 and load_connectivity_flag:
                # For efficiency reasons, connectivity is loaded only for the first mesh and copied (as a reference) to the rest.
                mesh = Cardiac3DMesh(vtk_path, load_connectivity_flag=True)
                self.meshes.append(mesh)
            else:
                mesh = Cardiac3DMesh(vtk_path, load_connectivity_flag=False)
                self.meshes.append(mesh)
                if load_connectivity_flag:
                    self.meshes[i].triangles = self.meshes[0].triangles

        if load_connectivity_flag:
            self.triangles = self.meshes[0].triangles

    def as_numpy_array(self):
        """ """
        kk = [x.points for x in self.meshes]
        try:
            kk = np.stack(kk, axis=0)
        except:
            # Handle this error better
            # embed()
            # raise ValueError(
            self.logger.error(
                """
            Not possible to create Numpy array for individual {id}. \
            The folder is likely to be incomplete.
            """.format(
                    id=self.subjectID
                )
            )
            raise ValueError

        return kk

    def _get_vtk_paths(self):

        # TODO: provide file pattern as argument to constructor
        # The current (hardcoded) file pattern is the one used by the SPASM output.
        fp = os.path.join(
            self._root_folder, "output/world2gimias/output.{time_frame}.vtk"
        )
        return [fp.format(time_frame=x) for x in self._time_frames_as_path]

    def _time_frames_to_int(self):
        """
        Convert cardiac phases like "ED" and "ES" (for end-diastole and end-systole) to the corresponding integer indices
        :return: None
        """
        self._time_frames_dict = {t: t for t in self.time_frames}
        self._time_frames_dict["ED"] = 1
        self._time_frames_dict["ES"] = self.ES_time_frame
        self.time_frames = [self._time_frames_dict[t] for t in self.time_frames]

    def __repr__(self):
        return "Time series of {} meshes (class {}) for subject {}".format(
            len(self.meshes), self.meshes[0].__class__.__name__, self.subjectID
        )

    def __getitem__(self, timeframe):
        return self.meshes[self._time_frames_dict[timeframe]]

    @property
    # TODO: TEST
    def ES_time_frame(self):
        with open(os.path.join(self._root_folder, "ES_time_step.csv"), "rt") as ff:
            self._ES_time_frame = int(ff.read().strip())
        return self._ES_time_frame

    @property
    # TODO: TEST
    def LVEF(self):
        with open(os.path.join(self._root_folder, "Ejection_fraction.csv"), "rt") as ff:
            self._LVEF = float(ff.read().strip())
        return self._LVEF

    @property
    # TODO: TEST
    def LVSV(self):
        with open(os.path.join(self._root_folder, "Stroke_volume.csv"), "rt") as ff:
            self._LVSV = float(ff.read().strip())
        return self._LVSV

    @property
    def _time_frames_as_path(self):
        return [
            "0" * (3 - len(str(t))) + str(t) for t in self.time_frames
        ]  # "001", "002", ..., "050"

    def generate_gif(self, gif_path, paraview_config):
        """
        Generate a GIF file showing the moving mesh, using Paraview.
        gif_path: path to the GIF output file
        paraview_config: object representing the Paraview config (specify what's needed)
        """
        raise NotImplementedError


# "~/data/PhD/meshes/vtk_meshes/2ch_full_cycle/1000215/output/world2gimias/output.001.vtk"


class CardiacMeshPopulation:

    """
    Class representing a population of cardiac meshes (either 3D or 4D),
    i.e. meshes for different individuals in a population

    Public attributes:
      meshes
      triangles
      subject_ids
      meanShape
      vertex_wise_stddev

    Usage example:
      mesh_pop = CardiacMeshPopulation(<ROOT_FOLDER>)
      mesh_pop[<SUBJECT_ID>] <--- either a Cardiac3DMesh or a Cardiac4DMesh object
    """

    def __init__(
        self,
        root_path=None,
        filename_pattern=None,
        time_frames=None,
        N_subj=None,
        shuffle=False,
        random_state=None,
        in_memory=True,
        logger=None,
    ):

        """
        #TODO: complete this docstring
        params:
            filename_pattern:
            time_frames:
            shuffle:
            random_state:
            N_subj:
            in_memory:
            logger:
        """

        self._root_path = root_path
        self._N_subj = N_subj
        self._shuffle = shuffle

        self._folders = [
            os.path.join(self._root_path, x) for x in os.listdir(self._root_path)
        ]

        if self._shuffle:
            random.seed(random_state)
            random.shuffle(self._folders)

        if self._N_subj is not None:
            self._folders = self._folders[: self._N_subj]

        self._logger = set_logger(logger)
        self.time_frames = time_frames

        # TODO: implement support for data accessing from disk directly
        if in_memory:
            self._load_data()
        else:
            raise NotImplementedError

    def _load_data(self):

        self.meshes, self.ids = [], []
        counter = 0

        for i, p in enumerate(tqdm(self._folders, unit="subjects")):
            try:
                c4dm = Cardiac4DMesh(p, time_frames=self.time_frames)

                try:
                    # TOFIX: this is very inefficient
                    # It's aimed to detect and bypass those cases where a VTK file is corrupt.
                    c4dm.as_numpy_array()
                except:
                    continue

                if i == 0:
                    self.time_frames = c4dm.time_frames

                id = c4dm.subjectID
                self.meshes.append(c4dm)
                self.ids.append(id)
                # counter += 1
                # if self._N_subj is not None and counter == self._N_subj:
                #    break
            except:
                # TODO: identify malformed folders
                self._logger.warning(
                    "Folder {} could not be read successfully".format(p)
                )

        # call triangles attribute from the Cardiac4DMesh class.
        self.triangles = self[0].triangles

    def __getitem__(self, indices):

        # TODO: implement an indexing scheme as the following
        """
        CMP = CardiacMeshPopulation(...)
        CMP[<SUBJECT_ID>]: a Cardiac4DMesh
        CMP[<SUBJECT_ID>, <TIMEFRAMES>] ---> CMP[<SUBJECT_ID>][<TIMEFRAMES>]
        CMP[<SUBJECT_ID>, <TIMEFRAMES>, [<PARTITIONS>]] ---> CMP[<SUBJECT_ID>,<TIMEFRAMES>][<PARTITIONS>]
        :param id:
        :return:
        """
        if isinstance(indices, int):
            int_idx = indices
            if int_idx >= 0 and int_idx < len(self.ids):
                return self.meshes[int_idx]
        elif isinstance(indices, str):
            # If a (single) string, it's interpreted as an individual's ID.
            subject_id = indices
            subject_index = self.ids.index(subject_id)
            return self.meshes[subject_index]
        elif isinstance(indices, tuple) or isinstance(indices, list):
            if len(indices) == 2:
                subject_id, timeframe = indices
                return self[subject_id][timeframe]
            elif len(indices) == 3:
                subject_id, timeframe, partition = indices
                return self[subject_id, timeframe][partition]

    @property
    def meanShape(self, mode=None):

        raise NotImplementedError
        return self._meanShape

    @property
    def vertex_wise_stddev(self, mode=None):
        raise NotImplementedError
        return self._stddev

    def _normalize(self):
        # Mean and std. are computed based on all the samples (not only the training ones). I think this makes sense.
        # Create self.is_normalized argument and set to True to track normalization status.
        self.mean, self.std = np.mean(self.point_clouds, axis=0), np.std(
            self.point_clouds, axis=0
        )
        self.point_clouds = (self.point_clouds - self.mean) / self.std
        self.is_normalized = True
        self._logger.info("Vertices normalized")

    def as_numpy_array(self):
        return np.stack([x.as_numpy_array() for x in self.meshes], axis=0)

    @property
    def shapePCA(self, n_comps=20):
        raise NotImplementedError
        try:
            self._shapePCA
        except AttributeError:
            # Code to implement shape PCA
            self._shapePCA = {"eigenvalues": eigenvals, "eigenvectors": eigenvecs}
            return self._shapePCA

    def generalisedProcrustes(self, scaling=True):

        if scaling:
            self._logger.info("Performing Procrustes analysis with scaling")
            self.reference_mesh = self.meshes[0]
            old_disparity, disparity = 0, 1  # random values
            it_count = 0
            while abs(old_disparity - disparity) / disparity > 1e-2 and disparity:
                old_disparity = disparity
                disparity = 0
                for i in range(len(self.point_clouds)):
                    # Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html
                    if self.procrustes_scaling:
                        mtx1, mtx2, _disparity = procrustes(
                            self.reference_mesh, self.point_clouds[i]
                        )
                        self.point_clouds[i] = np.array(
                            mtx2
                        )  # if self.procrustes_scaling else np.array(mtx1)
                    else:
                        # Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html
                        # Note that the arguments are swapped respect to the previous @procrustes function
                        R, s = orthogonal_procrustes(
                            self.point_clouds[i], self.reference_mesh
                        )
                        # Rotate
                        self.point_clouds[i] = np.dot(self.point_clouds[i], R)  # * s
                        # Mean point-wise MSE
                        _disparity = np.mean(
                            np.sqrt(
                                np.sum(
                                    np.square(
                                        self.point_clouds[i] - self.reference_mesh
                                    ),
                                    axis=1,
                                )
                            )
                        )
                    disparity += _disparity
                disparity /= self.point_clouds.shape[0]
                self.reference_mesh = self.point_clouds.mean(axis=0)
                it_count += 1
            self.procrustes_aligned = True
            self._logger.info(
                "Generalized Procrustes analysis with scaling performed after %s iterations"
                % it_count
            )

        else:
            self._logger.info("Performing Procrustes analysis without scaling")
            from scipy.linalg import orthogonal_procrustes

            self.reference_mesh = self.point_clouds[0]
            old_disparity, disparity = 0, 1

            # Center the meshes
            for i in range(len(self.point_clouds)):
                self.point_clouds[i] -= np.mean(self.point_clouds[i], 0)

                it_count = 0
                while abs(old_disparity - disparity) / disparity > 1e-4:
                    old_disparity = disparity
                    disparity = 0
                    for i in range(len(self.point_clouds)):
                        # Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html
                        R, s = orthogonal_procrustes(
                            self.point_clouds[i], self.reference_mesh
                        )

                        # Rotate
                        self.point_clouds[i] = np.dot(self.point_clouds[i], R)  # * s

                        # Mean point-wise MSE
                        _disparity = np.mean(
                            np.sqrt(
                                np.sum(
                                    np.square(
                                        self.point_clouds[i] - self.reference_mesh
                                    ),
                                    axis=1,
                                )
                            )
                        )

                        disparity += _disparity

                    disparity /= self.point_clouds.shape[0]
                    self.reference_mesh = self.point_clouds.mean(axis=0)
                    print(disparity)
                    it_count += 1

                self.procrustes_aligned = True
                self._logger.info(
                    "Generalized Procrustes analysis performed after %s iterations"
                    % it_count
                )

        raise NotImplementedError
        return rotation, translation
