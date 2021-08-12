import vtk
import numpy as np
import os
import meshio # tested with 2.3.0

'''
This module is aimed to simplify the implementation of common tasks on VTK triangular cardiac meshes,
that result overly convoluted if the usual VTK Python wrapper for C++ is used,
and render the code difficult to follow.
'''

LEGACY_2CHAMBER_SPASM = "LEGACY_2CHAMBER_SPASM"
LEGACY_4CHAMBER_MMF   = "LEGACY_4CHAMBER_MMF"
FULL_HEART_MODEL_MMF  = "FULL_HEART_MODEL_MMF"


class Cardiac3DMesh:

    '''
    This class represents a single cardiac mesh or point cloud.    
    '''
    
    def __init__(self, filename=None, subpartIDs=None, load_connectivity_flag=True, dataset_version=None):
        
        '''
        filename: path to a file for a VTK Polydata object.
        subpartIDs: either a label or a list of labels to subset the mesh for (e.g. "LV", "RV", etc.)
        load_connectivity_flag: boolean indicating whether or not to load the mesh connectivity information
        '''
        
        self.points = None        

        if filename is not None:
            self._filename = filename
            self._reader = vtk.vtkPolyDataReader()
            self._reader.SetFileName(self._filename)
            self._reader.Update()
                        
            self._load_point_cloud()                                    
            if load_connectivity_flag: self.load_connectivity()                
            self._load_partition_ids()
            self._infer_dataset_version()
        
        if subpartIDs is not None:
            self = self.extract_subpart(subpartIDs)

                        
    def _load_point_cloud(self):

        '''
        :return: numpy array where each element is a triple of (x, y, z)
        coordinates and a set of indices representing the links to that point
        '''

        output = self._reader.GetOutput()
        n_points = output.GetNumberOfPoints()        
        self.points = np.array([output.GetPoint(i) for i in range(n_points)])
        
           
    def load_connectivity(self, triangles=None, edges=None):

        '''
        triangles: if not None, must be a list of 3-tuples containing valid point indices
        edges: if not None, must be a list of 2-tuples containing valid point indices
        if both are None, the connectivity will be read from the VTK file
        '''

        if triangles is None and edges is None:            
            # TODO: raise an error if the VTK file has not been provided
            output = self._reader.GetOutput()
            self.n_cells = output.GetNumberOfCells()
            self.triangles = [[int(output.GetCell(l).GetPointId(k)) for k in (0, 1, 2)] for l in range(self.n_cells)]                
        elif triangles is not None:
            # TODO: perform better sanity check on the `triangles` argument
            if not isinstance(triangles, list):
                raise TypeError
            self.triangles = triangles
        else:
            self.triangles = self._triangles_from_edges(edges)            
            

    def _infer_dataset_version(self):

        if self.distinct_subparts == {1, 2, 4}:
            self._dataset_version = LEGACY_2CHAMBER_SPASM
        else:
            self._dataset_version = FULL_HEART_MODEL_MMF


    @property
    def edges(self):
        try:
            return self._edges
        except AttributeError:
            self._edges = self._edges_from_triangles(self.triangles)
            return self._edges

        
    def _edges_from_triangles(self, triangles):
        edges = []
        for x,y,z in triangles:            
            edges.extend([(x,y), (y,z), (x,z), (y,x), (z,y), (z,x)])
        return list(set(edges))
    
    
    def _triangles_from_edges(self):
        raise NotImplementedError
        #TODO: implement
        
        
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
        '''
        This works as an alias, for compatibility with third-party code.
        Probably it is better to use a wrapper (child) class for this.
        '''
        return self.points
    
    
    @property
    def f(self):
        '''
        This works as an alias, for compatibility with third-party code.
        Probably it is better to use a wrapper (child) class for this.
        '''
        return self.triangles

    
    def _load_partition_ids(self):

        '''
        Generate a list of the subpart IDs for each of the vertices (i.e. which partition of the mesh they belong to)
        '''

        output = self._reader.GetOutput()

        pp = output.GetPointData().GetAbstractArray(0)
        self.subpartID = [pp.GetValue(i) for i in range(self.n_points)]
        # self.subpartID = [int(pp.GetComponent(i, 0)) for i in range(self.n_points)]

        
    @property
    def n_points(self):
        '''
        Number of vertices in the mesh
        '''
        return len(self.points)
 

    @property
    def n_edges(self):
        '''
        Number of edges in the mesh (without double-counting)
        '''
        return len(self.edges)
    
    
    def extract_subpart(self, ids):
        '''
        ids: a label or a list of labels for the subpart/s to be extracted
        :return: a Cardiac3DMesh object representing the subpart to be extracted.
        '''
        
        ids = [ids] if not isinstance(ids, list) else ids

        subvtk = Cardiac3DMesh()
        subvtk.points = np.array([self.points[i] for i in range(self.n_points) if self.subpartID[i] in ids])
        subvtk.subpartID = np.array([self.subpartID[i] for i in range(self.n_points) if self.subpartID[i] in ids])

        if self.triangles is not None:
            point_ids = [ i for i, id in enumerate(self.subpartID) if id in ids ]
            point_ids_set = set(point_ids)
            triangles = [tuple(triangle) for triangle in self.triangles if all([pp in point_ids_set for pp in triangle])]

            id_mapping = { x:i for i, x in enumerate(point_ids) }
            subvtk.triangles = np.array([ tuple([id_mapping[x] for x in triangle]) for triangle in triangles ])

        return subvtk

      
    def __repr__(self):
        return "Point cloud\n\n {} \n\n with connectivity\n\n{}".format(self.points.__str__(), self.triangles.__str__())
        

    def __getitem__(self, id):
        
        if self._dataset_version == LEGACY_2CHAMBER_SPASM:
            if id == "LV_endo":
                return self.extract_subpart(1)
            elif id == "LV_epi":
                return self.extract_subpart(2)
            elif id == "LV":
                return self.extract_subpart([1,2])
            elif id == "RV_endo" or id == "RV":
                return self.extract_subpart(4)
            else:
                raise ValueError("{} is not a valid partition (use LV_endo, LV_epi, LV, RV or RV_endo).".format(id))
        elif self._dataset_version == FULL_HEART_MODEL_MMF:
            id = [id] if not isinstance(id, list) else id
            if all([x in self.distinct_subparts for x in id]):
                return self.extract_subpart(id)
            else:
                raise ValueError("{} is not a valid partition (use {} or combinations thereof)".format(id, ", ".join(sorted(list(self.distinct_subparts)))))

            
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
            self._adj_matrix = sp.csc_matrix((
            np.ones(len(self.edges)), 
              ([x[0] for x in self.edges], 
               [x[1] for x in self.edges])
            ))
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
            cells={'triangle': np.array(self.triangles)},
            point_data={'subpartID': self.subpartID}
        )                
