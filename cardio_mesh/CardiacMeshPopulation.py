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
        self._random_state = random_state

        self._folders = [
            os.path.join(self._root_path, x) for x in os.listdir(self._root_path)
        ]

        if self._shuffle:
            random.seed(self._random_state)
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