import vtk
import numpy as np
import os
import meshio # tested with 2.3.0

'''
This module is aimed to simplify the implementation of common tasks on VTK meshes,
that result overly convoluted if the usual VTK Python wrapper for C++ is used,
and render the code difficult to follow.
'''

class VTKObject():

    def __init__(self, filename=None, subpartIDs=None, load_connectivity=True):

        self.filename = filename
        self.n_points = None
        self.n_cells = None
        self.points = None
        self.edges = []
        self.triangles = []
        self.neighbors_dict = {}
        self.load_connectivity = load_connectivity

        if filename is not None:
            self.reader = vtk.vtkPolyDataReader()
            self.reader.SetFileName(self.filename)
            self.reader.Update()
            self.getMesh()
            # self.extractTriangles()
            self.extractPartitionIDs()
            self.v = self.points # for compatibility with other code
            self.f = np.array(self.triangles) # idem
        
        if subpartIDs is not None:
            self = self.extractSubpart(subpartIDs)


    def getMesh(self):

        '''
        :return: numpy array where each element is a triple of (x, y, z)
        coordinates and a set of indices representing the links to that point
        '''

        output = self.reader.GetOutput()
        self.n_points = output.GetNumberOfPoints()
        self.n_cells = output.GetNumberOfCells()
        self.points = np.array([output.GetPoint(i) for i in range(self.n_points)])
        if self.load_connectivity:
            self.loadConnectivity()
       
    
    def loadConnectivity(self):
 
        output = self.reader.GetOutput()
        for i in range(self.n_cells):
            pts_cell = output.GetCell(i)
            for j in range(pts_cell.GetNumberOfEdges()):
                self.edges.append(([int(pts_cell.GetEdge(j).GetPointId(i)) for i in (0, 1)]))
                self.neighbors_dict.get(self.edges[-1][0], []).append(self.edges[-1][1])
        self.triangles = [[int(output.GetCell(j).GetPointId(i)) for i in (0, 1, 2)] for j in range(self.n_cells)]


    def extractTriangles(self):
        output = self.reader.GetOutput()
        self.triangles = [[int(output.GetCell(j).GetPointId(i)) for i in (0, 1, 2)] for j in range(self.n_cells)]


    def extractPartitionIDs(self):

        '''
        Generate a list of the subpart IDs for each of the vertices (i.e. which partition of the mesh they belong to)
        '''

        output = self.reader.GetOutput()

        pp = output.GetPointData().GetAbstractArray(0)
        self.subpartID = [pp.GetValue(i) for i in range(self.n_points)]


    def extractSubpart(self, ids):

        subvtk = VTKObject()
        subvtk.points = np.array([self.points[i] for i, x in enumerate(self.points) if self.subpartID[i] in ids])
        subvtk.subpartID = np.array([self.subpartID[i] for i, x in enumerate(self.points) if self.subpartID[i] in ids])
        subvtk.n_points = len(subvtk.points)

        if self.edges is not None:
            point_ids = [ i for i, id in enumerate(self.subpartID) if id in ids ]
            point_ids_set = set(point_ids)
            triangles = [tuple(triang) for triang in self.triangles if all([pp in point_ids_set for pp in triang])]

            id_mapping = { x:i for i, x in enumerate(point_ids) }
            subvtk.triangles = np.array([ tuple([id_mapping[x] for x in tr]) for tr in triangles ])

            subvtk.v = subvtk.points
            subvtk.f = subvtk.triangles

        return subvtk

    def get_adj_matrix(self):

        """
        Returns a sparse matrix (of size #verts x #verts) where each nonzero
        element indicates a neighborhood relation. For example, if there is a
        nonzero element in position (15,12), that means vertex 15 is connected
        by an edge to vertex 12.
        """

        from scipy import sparse as sp

        adj_matrix = sp.csc_matrix((
            np.ones(len(self.edges)), 
              ([x[0] for x in self.edges], 
               [x[1] for x in self.edges])
            ))

        return adj_matrix

    def adj_matrix_to_edges(self, adj_matrix):
        from scipy import sparse as sp
        non_zero_indices = sp.find(adj_matrix)
        return zip(non_zero_indices[0], non_zero_indices[1])
    

    # mesh to vtk
    # TODO: Add subpartID as a label in the output file
    # meshio doesn't seem to support strings as point data
    def SaveMeshToVTK(self, filename):
        #meshio.write_points_cells(
        #    filename,
        #    self.points,
        #    cells={'triangle': np.array(self.triangles)}
        #)
        mesh = meshio.Mesh(self.v, {'triangle':self.f})
        meshio.write(filename, mesh)

