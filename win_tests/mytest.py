import os
print("Hi")
import pyxyz as m
print("Hi!!!")
import pyxyz.test
pyxyz.test.run_tests()


if __name__ == "__main__":
    x = m.Confpool()
    print("Hi")
    x.include_from_file("aminoacid_ensemble.xyz", energy=lambda x: float(x.strip()))
    # x.energy_filter(1.0, etype='kcal/mol')
    # x.include("testcs.xyz")
    # x.update_description(lambda e, d: "Energy = " + repr(e))
    x.distance_filter(19, 1, lambda x: x > 1.7)
    # x.valence_filter(2, 1, 31, lambda x: x > 160.0)
    # x.dihedral_filter(7,9,10,11, lambda x: abs(x - 116.0) < 10.0)
    x.sort()
    x.save('test_res2.xyz')
