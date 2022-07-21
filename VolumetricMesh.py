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
    This class represents a single cardiac mesh or point cloud, for volumetric meshes without subparts
    """

    def __init__(
        self,
        filename=None,
        load_connectivity_flag=True,
        logger=None,
    ):

        """
        filename: path to a file for a VTK Polydata object or to a PKL.
        load_connectivity_flag: boolean indicating whether or not to load the mesh connectivity information
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

            elif self._filename.endswith(".pkl"):

                with open(self._filename, "rb") as f:
                    dict = pkl.load(f)
                
                # We can assume vertices and faces are numpy arrays
                self.points = dict['points']
                self.triangles = dict['triangles']


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
        dict = {"points" : self.points, "triangles" : self.triangles}

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