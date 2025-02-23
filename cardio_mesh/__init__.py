# For relative imports to work

import os, sys; 

sys.path.append(BASE_DIR := os.path.dirname(os.path.realpath(__file__)))

# import CardiacMesh
from CardiacMesh import Cardiac3DMesh

import procrustes
