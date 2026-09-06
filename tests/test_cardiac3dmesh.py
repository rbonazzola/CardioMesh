import os

import numpy as np
import pytest

from CardiacMesh import Cardiac3DMesh

SAMPLE_DATA = os.path.join(os.path.dirname(__file__), "sample_data")
FHM_VTK = os.path.join(SAMPLE_DATA, "full_heart_model", "full_heart_model.vtk")
LEGACY_VTK = os.path.join(SAMPLE_DATA, "legacy_biventricular_model", "output.001.vtk")


@pytest.fixture(scope="module")
def fhm_mesh():
    return Cardiac3DMesh(FHM_VTK)


@pytest.fixture(scope="module")
def legacy_mesh():
    return Cardiac3DMesh(LEGACY_VTK)


class TestVtkLoader:
    def test_missing_file_raises(self):
        with pytest.raises(FileExistsError):
            Cardiac3DMesh("does_not_exist.vtk")

    def test_full_heart_model_loads(self, fhm_mesh):
        assert fhm_mesh._dataset_version == "FULL_HEART_MODEL_MMF"
        assert fhm_mesh.n_points > 0
        assert fhm_mesh.triangles.shape[1] == 3
        assert len(fhm_mesh.subpartID) == fhm_mesh.n_points
        assert fhm_mesh.distinct_subparts >= {"LV", "RV", "LA", "RA"}

    def test_legacy_biventricular_model_loads(self, legacy_mesh):
        assert legacy_mesh._dataset_version == "LEGACY_2CHAMBER_SPASM"
        assert legacy_mesh.distinct_subparts == {1, 2, 4}

    def test_aliases_and_len(self, legacy_mesh):
        assert legacy_mesh.v is legacy_mesh.points
        assert legacy_mesh.f is legacy_mesh.triangles
        assert len(legacy_mesh) == legacy_mesh.n_points


class TestSubpartExtraction:
    def test_extract_single_partition(self, fhm_mesh):
        lv = fhm_mesh["LV"]
        assert 0 < lv.n_points < fhm_mesh.n_points
        assert set(lv.subpartID) == {"LV"}
        # every triangle must reference valid (remapped) point indices
        assert lv.triangles.max() < lv.n_points

    def test_extract_multiple_partitions(self, fhm_mesh):
        combined = fhm_mesh["LV", "RV"]
        assert set(combined.subpartID) == {"LV", "RV"}
        assert combined.n_points == fhm_mesh["LV"].n_points + fhm_mesh["RV"].n_points

    def test_extract_from_legacy_mesh(self, legacy_mesh):
        # Regression test: extracting a subpart used to crash with
        # AttributeError because _subpart_id_mapping_str_to_int/_int_to_str
        # were only set for the FULL_HEART_MODEL_MMF dataset version.
        lv = legacy_mesh["LV"]
        assert 0 < lv.n_points < legacy_mesh.n_points
        assert set(lv.subpartID) == {1, 2}

    def test_invalid_partition_raises(self, fhm_mesh):
        with pytest.raises(ValueError):
            fhm_mesh["not_a_real_partition"]


class TestTopology:
    def test_edges_are_symmetric_and_deduplicated(self, legacy_mesh):
        lv = legacy_mesh["LV"]
        edges = lv.edges
        assert lv.n_edges == len(edges)
        edge_set = set(edges)
        # every edge (a, b) must have its reverse (b, a) present
        assert all((b, a) in edge_set for a, b in edges)

    def test_neighbors_dict(self, legacy_mesh):
        lv = legacy_mesh["LV"]
        neighbors = lv.neighbors_dict
        assert len(neighbors) > 0
        # each recorded neighbor must correspond to a real edge
        a, b = lv.edges[0]
        assert b in neighbors[a]

    def test_adj_matrix_matches_edge_count(self, legacy_mesh):
        lv = legacy_mesh["LV"]
        assert lv.adj_matrix.nnz == lv.n_edges


class TestPklRoundTrip:
    def test_save_and_reload(self, legacy_mesh, tmp_path):
        lv = legacy_mesh["LV"]
        pkl_path = tmp_path / "lv.pkl"
        lv.save_to_pkl(str(pkl_path))

        reloaded = Cardiac3DMesh(str(pkl_path))

        assert np.array_equal(reloaded.points, lv.points)
        assert np.array_equal(reloaded.triangles, lv.triangles)
        assert np.array_equal(reloaded.subpartID, lv.subpartID)


class TestNpyLoader:
    def _write_triangle_fixture(self, tmp_path):
        # A minimal two-triangle quad, split into "A" and "B" subparts.
        points = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        )
        triangles = [(0, 1, 2), (0, 2, 3)]
        subpart_ids = ["A", "A", "B", "B"]

        points_path = tmp_path / "points.npy"
        faces_path = tmp_path / "faces.csv"
        subpart_path = tmp_path / "subpart.txt"

        np.save(points_path, points)
        with open(faces_path, "w") as f:
            f.write("\n".join(",".join(str(i) for i in t) for t in triangles))
        with open(subpart_path, "w") as f:
            f.write("\n".join(subpart_ids))

        return points, triangles, points_path, faces_path, subpart_path

    def test_loads_points_faces_and_subpart_ids(self, tmp_path):
        points, triangles, points_path, faces_path, subpart_path = (
            self._write_triangle_fixture(tmp_path)
        )

        mesh = Cardiac3DMesh(
            str(points_path),
            faces_filename=str(faces_path),
            subpart_id_filename=str(subpart_path),
        )

        assert np.array_equal(mesh.points, points)
        assert np.array_equal(mesh.triangles, np.array(triangles))
        assert mesh.subpartID == ["A", "A", "B", "B"]

    def test_subpart_id_count_mismatch_raises(self, tmp_path):
        _, _, points_path, faces_path, subpart_path = self._write_triangle_fixture(
            tmp_path
        )
        # Truncate the subpart ID file so it no longer matches the point count.
        subpart_path.write_text("A\nA\n")

        with pytest.raises(ValueError):
            Cardiac3DMesh(
                str(points_path),
                faces_filename=str(faces_path),
                subpart_id_filename=str(subpart_path),
            )
