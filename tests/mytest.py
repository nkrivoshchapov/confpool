import os
import confpool as m
# import pyxyz.test
# pyxyz.test.run_tests()

def func(x):
    # print("x is " + repr(x))
    # print("returns " + repr(x["OH contact"]))
    return x["OH contact"] < 2.0

if __name__ == "__main__":
    x = m.Confpool()

    curdir = os.path.dirname(os.path.realpath(__file__))
    x.include_from_file(os.path.join(curdir, "hydrogen_atom.xyz"))
    x.key_from_description("Energy", lambda x: float(x.strip().split('=')[1]))
    print("Size = {}".format(x.size()))

    x = m.Confpool()
    x.include_from_file(os.path.join(curdir, "aminoacid_ensemble.xyz"))
    x.key_from_description("Energy", lambda x: float(x.strip()))
    print("XXXSize = {}".format(x.size()))
    x.distance_to_key("OH contact", 1, 19)
    x.vangle_to_key("VA", 7, 19, 1)
    x.dihedral_to_key("TA", 8, 6, 7, 19)
    # x.filter(lambda x: abs(x['TA']) > 160.0)
    # x.sort("VA")
    # x.filter(lambda x: x["OH contact"] < 2.0)
    print("A Count = {}".format(x.count(func)))
    print("B Count = {}".format(x.count(lambda mol: mol.l(1, 19) < 2.0)))
    # x.filter(func)
    # x.upper_cutoff("OH contact", 2.0)
    # x.lower_cutoff("OH contact", 2.0)
    # x.update_description(lambda x: "D = {}. Positive energy = {} (old = {})".format(x["OH contact"], -x["Energy"], x.descr().strip()))
    x.update_description(lambda x: "TA = {} VA = {}".format(x["TA"], x["VA"]))
    x.rmsd_filter(0.3)
    print(repr(x[0].xyz()))
    x.save('lolkek.xyz')
    print("Count = {}".format(x.count(func)))
    print("XXXSize = {}".format(x.size()))

    raise Exception("pause")
    # x.energy_filter(1.0, etype='kcal/mol')
    # x.include_from_file("testcs.xyz")
    x.update_description(lambda e, d: "Energy = " + repr(e))
    # x.distance_filter(3, 31, lambda x: x > 1.7)
    # x.valence_filter(2, 1, 31, lambda x: x > 160.0)
    x.dihedral_filter(7,9,10,11, lambda x: abs(x - 116.0) < 10.0)
    x.sort()
    print("Size = " + repr(x.size()))
    # for i in range(x.size()):
    #     print("i = {} -- {}".format(i, repr(x.get_structure(i))))
    print("Syms = " + repr(x.get_atom_symbols()))
    x.save('test_res2.xyz')
