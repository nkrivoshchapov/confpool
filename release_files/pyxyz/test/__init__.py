import os
from .. import Confpool

def run_tests():
    # from pyxyz import Confpool
    H2KC = 627.509474063
    KC2H = 1 / H2KC
    curdir = os.path.dirname(os.path.realpath(__file__))
    print("Running tests...")

    p = Confpool()
    p.include_from_file(os.path.join(curdir, "hydrogen_atom.xyz"))
    assert p.size == 1

    p = Confpool()
    p.include_from_file(os.path.join(curdir, "aminoacid_single.xyz"))
    p.filter(lambda m: m.l(1, 19) < 2.0)
    assert p.size == 0

    p = Confpool()
    p.include_from_file(os.path.join(curdir, "aminoacid_ensemble.xyz"))
    p["Energy"] = lambda m: float(m.descr.strip())
    p.upper_cutoff("Energy", 5.0 * KC2H)
    p.filter(lambda m: m.l(1, 19) < 2.0)
    p.sort("Energy")
    p.descr = lambda m: "Energy = {}".format(m["Energy"])
    p.save(os.path.join(curdir, "aminoacid_ensemble_current.xyz"))
    assert p.size == 2
    assert len(repr(p[0].xyz)) > 0
    assert len(repr(p[1].xyz)) > 0
    ref_data = open(os.path.join(curdir, "aminoacid_ensemble_ref.xyz"), 'r').read()
    current_data = open(os.path.join(curdir, "aminoacid_ensemble_ref.xyz"), 'r').read()
    assert ref_data == current_data

    p = Confpool()
    p.include_from_file(os.path.join(curdir, "crest_conformersA.xyz"))
    p.include_from_file(os.path.join(curdir, "crest_conformersB.xyz"))
    assert p.size == 365
    n_del = 0
    p["Energy"] = lambda m: float(m.descr.strip())
    n_del += p.upper_cutoff("Energy", 5.0 * KC2H)
    assert n_del == 62
    p.sort("Energy")
    n_del += p.rmsd_filter(0.3)
    p.descr = lambda m: "Conf #{}; Energy = {:6f} a.u.".format(m.idx + 1, m["Energy"])
    assert n_del == 333

    p["H2O dist"] = lambda m: m.l(26, 27)
    for i in range(p.size):
        assert isinstance(p[i].l(26, 27), float)
    assert p.count(lambda m: m["H2O dist"] < 1.8) == 25

    def quality_check(m):
        if m["H2O dist"] < 1.8:
            comment = "Nice"
        else:
            comment = "Garbage"
        return "{}; E={:3f} a.u.".format(comment, m["Energy"])
    p.descr = quality_check
    assert p[0].descr == "Nice; E=-45.344433 a.u."

    hbond_condition = lambda m: m.l(25, 10) < 2.0 and m.v(24, 25, 10) > 150.0 and abs(m.z(26, 23, 24, 25)) > 120.0
    num = p.count(hbond_condition)
    assert num == 5
    p.filter(hbond_condition)

    min_ener = min(p["Energy"]) # p["Energy"] creates a list of floats
    p["E.Rel"] = lambda m: H2KC * (m["Energy"] - min_ener)
    p.descr = lambda m: "{}; E(rel.)={:3f} kcal/mol".format(m.descr, m["E.Rel"])
    p.save('check.xyz')

    p = Confpool()
    p.include_from_file(os.path.join(curdir, "crest_conformersA.xyz"))
    p.include_from_file(os.path.join(curdir, "crest_conformersB.xyz"))
    p["Energy"] = lambda m: float(m.descr.strip())
    p.sort("Energy")
    p.rmsd_filter(0.3)
    for i in range(p.size):
        if hbond_condition(p[i]) and p[i].l(26, 27) < 1.8:
            p[i]["MyCondition"] = 1.0
        else:
            p[i]["MyCondition"] = 0.0
    assert isinstance(p.as_table(), dict)

    print("Tests were completed successfully.")