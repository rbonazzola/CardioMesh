LEGACY_2CHAMBER_SPASM = "LEGACY_2CHAMBER_SPASM"

LEGACY_2CH_SUBPART_IDS = {     
    # "all": None,
    "LV": [1, 2],
    "LV_endo": [1],
    "LV_epi": [2],
    "RV_endo": [4],
    "RV": [4]
}

LEGACY_4CHAMBER_MMF   = "LEGACY_4CHAMBER_MMF"

LEGACY_4CH_SUBPART_IDS = {     
    # "all": None,
    "LV": [1, 2],
    "LV_endo": [1],
    "LV_epi": [2],
    "LA": [3],
    "RV_endo": [4],
    "RV": [4],
    "RA": [5]
}

FULL_HEART_MODEL_MMF  = "FULL_HEART_MODEL_MMF"
FHM_SUBPART_AS_INT = {
    'LV': 1,
    'RV': 2,
    'RA': 3,
    'LA': 4,
    'aorta': 5,
    'PA': 6,
    'TVP': 7,
    'AVP': 8,
    'MVP': 9,
    'PVP': 10,
    'PV1': 11,
    'PV7': 12,
    'PV6': 13,
    'PV4': 14,
    'PV5': 15,
    'PV2': 16,
    'PV3': 17,
}

FHM_SUBPART_AS_STR = {v: k for k,v in FHM_SUBPART_AS_INT.items()}

FHM_SUBPART_IDS = {
    "AO":["aorta"],
}

closed_partitions = {
  "LA_closed" : ("LA", "MVP", "PV1", "PV2", "PV3", "PV4", "PV5"),
  "RA_cloase" : ("RA", "TVP", "PV6", "PV7"),
  "LV_closed" : ("LV", "AVP", "MVP"),
  "RV_closed" : ("RV", "PVP", "TVP"),
  "BV_closed" : ("LV", "AVP", "MVP", "RV", "PVP", "TVP"),
  "aorta" : ("aorta",)
}
