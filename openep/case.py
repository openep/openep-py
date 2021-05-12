import numpy as np
from typing import Any, Dict


# # Print the current working directory
# print("Current working directory: {0}".format(cwd))



class Case:
    def __init__(
        self,
        name: str,
        nodes: np.array,
        indices: np.array,
        fields: Dict[str, np.array],
        electric: Dict[str, np.array],
        rf: Dict[str, np.array],
        other_data: Dict[str, Any] = {},
    ):
        self.name: str = name
        self.nodes: np.array = nodes
        self.indices: np.array = indices
        self.fields: Dict[str, np.array] = dict(fields)
        self.electric: Dict[str, np.array] = dict(electric)
        self.rf: Dict[str, np.array] = dict(electric)
        self.other_data: Dict[str, Any] = dict(other_data)

    def __repr__(self):
        return f"{self.name}( nodes: {self.nodes.shape} indices: {self.indices.shape} fields: {tuple(self.fields)} )"

    def create_mesh(self, vertex_norms=True, recenter=True):
        inds_inverted = self.indices[:, [0, 2, 1]]
        inds = np.vstack([self.indices, inds_inverted])  # include each triangle twice for both surfaces

        mesh = trimesh.Trimesh(self.nodes, inds, process=False)
        mesh._kwargs["parent"]=self

        if vertex_norms:
            mesh.vertex_normals  # compute vertex normals

        if recenter:
            mesh.apply_translation(-mesh.centroid)  # recenter mesh to origin, helps with lighting in default scene

        return mesh


def load_case(filename, name=None, exclude_pattern=".*#.*"):
    """
    Load a Case object from the given file. This assumes a number of names for objects found in the file.
    """
    dat = load_mat(filename, exclude_pattern)

    if name is None:
        name = os.path.basename(filename)

    nodes = dat["userdata/surface/triRep/X"].T
    inds = dat["userdata/surface/triRep/Triangulation"].T - 1
    act, bip = dat["userdata/surface/act_bip"]
    uni, imp, frc = dat["userdata/surface/uni_imp_frc"]

    fields = dict(act=act, bip=bip, uni=uni, imp=imp, frc=frc)
    electric = {}
    rf = {}
    other_data = {}

    for k, v in dat.items():
        _, section, *objname = k.split("/", 2)

        if objname and section in ("electric", "rf"):
            dobj = rf if section == "rf" else electric
            dobj[objname[0]] = v
        else:
            other_data[k] = v

    return Case(name, nodes, inds, fields, electric, rf, other_data)


case = load_case("new_dataset_1.mat")
print(case, tuple(case.rf))
