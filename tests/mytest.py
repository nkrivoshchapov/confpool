import pyxyz as m
import pyxyz.test
pyxyz.test.run_tests()

if __name__ == "__main__":
    # x = m.Confpool()
    # x.include_from_file("aminoacid_ensemble_ref.xyz", energy=lambda x: float(x.strip().split('=')[1]))
    # x.save('checkout.xyz')
    x = m.Confpool()
    x.include_from_file("testcs.xyz", energy=lambda x: float(x.strip()))
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
