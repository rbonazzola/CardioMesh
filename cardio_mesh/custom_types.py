from typing import Dict, List, Optional, Hashable, Literal, Union
import numpy as np

PointCloud = np.ndarray
DynamicPointCloud = np.ndarray

PointCloudID = Hashable

Transformations = Dict[Literal["traslation", "rotation"], np.ndarray]
ProcrustesResult = Dict[PointCloudID, Transformations]

PartitionID = Literal["LV", "RV", "LA", "RA", "LV_closed", "RV_closed", "LA_closed", "RA_closed"]