import os
from .. import Confpool

H2KC = 627.509474063

def run_tests():
    print("Running tests...")
    curdir = os.path.dirname(os.path.realpath(__file__))
    x = Confpool()
    x.include_from_file(os.path.join(curdir, "hydrogen_atom.xyz"))
    x.key_from_description("Energy", lambda x: float(x.strip().split('=')[1]))
    assert x.size() == 1

    x = Confpool()
    x.include_from_file(os.path.join(curdir, "aminoacid_single.xyz"))
    x.distance_to_key("OH contact", 1, 19)
    x.filter(lambda x: x["OH contact"] < 2.0)
    assert x.size() == 0

    x = Confpool()
    x.include_from_file(os.path.join(curdir, "aminoacid_ensemble.xyz"))
    x.key_from_description("Energy", lambda x: float(x.strip()))
    x.update_description(lambda x: "Energy = " + repr(x['Energy']))
    x.filter(lambda x: x["Energy"] < 5.0 / H2KC)
    x.distance_to_key("OH contact", 1, 19)
    x.filter(lambda x: x["OH contact"] < 2.0)
    x.sort("Energy")
    x.save(os.path.join(curdir, "aminoacid_ensemble_current.xyz"))
    assert x.size() == 2
    assert len(repr(x.get_structure(0))) > 0
    assert len(repr(x.get_structure(1))) > 0
    ref_data = open(os.path.join(curdir, "aminoacid_ensemble_ref.xyz"), 'r').read()
    current_data = open(os.path.join(curdir, "aminoacid_ensemble_ref.xyz"), 'r').read()
    assert ref_data == current_data
    print("Tests were completed successfully.")
