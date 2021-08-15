import vtk
import numpy as np
import os
import meshio # tested with 2.3.0
from Constants import *

'''
This module is aimed to simplify the implementation of common tasks on VTK triangular meshes,
that result overly convoluted if the usual VTK Python wrapper for C++ is used,
and render the code difficult to follow.
'''

class Cardiac3DMesh:

    '''
    This class represents a single cardiac mesh or point cloud.    
    '''
    
    def __init__(self, filename=None, subpartIDs=None, load_connectivity_flag=True, dataset_version=None):
        
        '''
        filename: path to a file for a VTK Polydata object.
        subpartIDs: either a label or a list of labels to subset the mesh for (e.g. "LV", "RV", etc.)
        load_connectivity_flag: boolean indicating whether or not to load the mesh connectivity information
        dataset_version: if None, it is inferred automatically based on the subpart IDs. Currently supported values are "LEGACY_2CHAMBER_SPASM" and "FULL_HEART_MODEL_MMF"
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

        '''
        Infer the dataset version. Currently, "LEGACY_2CHAMBER_SPASM" and "FULL_HEART_MODEL_MMF" are supported.
        '''

        #TODO: decide if there is a better to infer the partition
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
    

    def __repr__(self):
        return "Point cloud\n\n {} \n\n with connectivity\n\n{}".format(self.points.__str__(), self.triangles.__str__())
        

    def show(self):
        from trimesh import Trimesh
        return Trimesh(self.v, self.f).show()
        

    def _extract_subpart(self, ids):
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


    def _map_subpart_ids(self, ids):
        '''

        '''

        ids = list(ids) if isinstance(ids, tuple) else ids
        ids = [ids] if not isinstance(ids, list) else ids

        # 1 -> [1]
        # 1,2 -> [1,2]
        # "LV" -> ["LV"]
        # "LV", "RV"  -> ["LV", "RV"]
        
        possible_values = list(self.distinct_subparts)
        possible_values += [ x for x in self._subpart_id_mapping if all([y in self.distinct_subparts for y in self._subpart_id_mapping[x]])] 

        kk = []

        for id in ids:
            for x in self._subpart_id_mapping.get(id, [id]):
                if x in possible_values:                
                    kk.append(x)
                else:
                    raise ValueError("{} is not a valid partition (use {} or combinations thereof)".format(x, ", ".join(sorted([str(x) for x in possible_values]))))

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


class Cardiac4DMesh:

    '''
    Class representing a collection of cardiac meshes for one individual, across the cardiac cycle.
    '''

    def __init__(self, root_folder, time_frames=None):
        
        '''
        root_folder: path to the PyCardioX output for the given individual      
        time_frames: list|tuple containing 1-50 and/or "ED"/"ES"
        '''
        
        self._root_folder = root_folder        
        if time_frames is None:            
            self.time_frames = [i+1 for i in range(50)]
        else:
            self.time_frames = time_frames
        self._time_frames_to_int()                
        self._vtk_paths = self._get_vtk_paths()                
        self._load_meshes()
        
        
    @property
    def subjectID(self):
        self._subjectID = os.path.basename(self._root_folder.strip("/"))
        return self._subjectID
        
    
    def _load_meshes(self, load_connectivity_flag=True):
        self.meshes = []
        for i, vtk_path in enumerate(self._vtk_paths):
             if i == 0 and load_connectivity_flag:
                 self.meshes.append(Cardiac3DMesh(vtk_path, load_connectivity_flag=True))
             else:
                 self.meshes.append(Cardiac3DMesh(vtk_path, load_connectivity_flag=False))
                 if load_connectivity_flag:
                     self.meshes[i].triangles = self.meshes[0].triangles
        
        if load_connectivity_flag:
            self.triangles = self.meshes[0].triangles
    
    
    def _get_vtk_paths(self):
        fp = os.path.join(self._root_folder, "output/world2gimias/output.{time_frame}.vtk")
        return [ fp.format(time_frame=x) for x in self._time_frames_as_path]
        
    
    def _time_frames_to_int(self):
        self._time_frames_dict = {t:t for t in self.time_frames}
        self._time_frames_dict["ED"] = 1
        self._time_frames_dict["ES"] = self.ES_time_frame        
        self.time_frames = [self._time_frames_dict[t] for t in self.time_frames]

        
    def __repr__(self):
        return "Time series of {} meshes (class {}) for subject {}.".format(len(self.meshes), self.meshes[0].__class__.__name__, self.subjectID)

    
    def __getitem__(self, timeframe):
        return self.meshes[timeframe]            
          
                                    
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
        return ["0"*(3-len(str(t))) + str(t) for t in self.time_frames] # "001", "002", ..., "050"
    
    
    def generate_gif(self, gif_path, paraview_config):
        '''
        Generate a GIF file showing the moving mesh, using Paraview.
        gif_path: path to the GIF output file
        paraview_config: object representing the Paraview config (specify what's needed)
        '''
        raise NotImplementedError
        #TODO: implement
                        

class CardiacMeshPopulation:
    
    def __init__(self, root_path, time_frames=("ED", "ES")):
        
        self._root_path = root_path
        self.CardiacPopulation = None
        self.subjectIDs = None
        raise NotImplementedError

        
    def __getitem__(self, id):     
        raise NotImplementedError
        idx = self.subjectIDs.index(id)
        return self.CardiacPopulation[idx]
    
    
    def generalisedProcrustes(self, scaling=True):       
        raise NotImplementedError
        return rotation, translation
    
    
    @property
    def meanShape(self, mode=None):        
        raise NotImplementedError        
        return self._meanShape
    
    
    @property
    def vertex_wise_stddev(self, mode=None):        
        raise NotImplementedError        
        return self._stddev
    
    
    @property    
    def shapePCA(self, n_comps=20):        
        raise NotImplementedError
        try:
            self._shapePCA
        except AttributeError:
            # Code to implement shape PCA
            self._shapePCA = {"eigenvalues":eigenvals, "eigenvectors": eigenvecs}
            return self._shapePCA
