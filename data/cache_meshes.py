import os, sys
import shlex
from subprocess import check_output

# Move to the root folder of the Git repository
os.chdir(check_output(shlex.split("git rev-parse --show-toplevel")).strip().decode('ascii'))
sys.path.append(f"{os.environ['HOME']}/01_repos/CardioMesh")

CARDIAC_COMA_REPO = f"{os.environ['HOME']}/01_repos/CardiacCOMA"

import pandas as pd
import pickle as pkl
import numpy as np
from pprint import pprint

import CardiacMesh
from CardiacMesh import Cardiac3DMesh

end_systole_timeframes = pd.read_csv(f"{CARDIAC_COMA_REPO}/data/cardio/end_systole_timeframes.csv")
end_systole_timeframes = dict(zip(end_systole_timeframes.id, end_systole_timeframes.end_systole_index))

root_folder = f"{CARDIAC_COMA_REPO}/data/cardio/meshes"

N = 40000
ids = [int(id) for id in os.listdir(root_folder)]
ids = ids[:N]


FACES = "data/faces_fhm_10pct_decimation.csv"
SUBPART_IDS = "data/subpartIDs_FHM_10pct.txt"
PROCRUSTES_FILE = "data/procrustes_transforms_FHM_35k.pkl"
procrustes_transforms = pkl.load(open(PROCRUSTES_FILE, "rb"))

os.makedirs("data/cached", exist_ok=True)


for part in ["LV", "RV", ("LV", "RV")]:
    
    for phase in ["ED", "ES"]:
                                
        meshes = {}
        
        for id in ids:
        
            # phases: ED and ES.            
            if phase == "ED":
                timeframe = "001"
            elif phase == "ES":
                timeframe = str(end_systole_timeframes[id]).zfill(3)  
            
            npy_file = f"{root_folder}/{id}/models/FHM_time{timeframe}.npy"            
            
            try:
                mesh = Cardiac3DMesh(
                  filename=npy_file,
                  faces_filename=FACES,
                  subpart_id_filename=SUBPART_IDS                                
                )
            except FileNotFoundError:
                pass
            
            # pc  = np.load(npy_file)
                  
            mesh = mesh[part]
            mesh.points = CardiacMesh.transform_mesh(mesh.points, **procrustes_transforms[str(id)])
            meshes[id] = mesh.points

        # convert ("LV", "RV") into LVRV
        file_prefix = f"{''.join(part) if isinstance(part, tuple) else part}{phase}"
        pkl_file = f"data/cached/{file_prefix}_{len(meshes)}subjects.pkl"
        pkl.dump(meshes, open(pkl_file, "wb"))