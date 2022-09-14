import os
from pyxyz import Confpool
import numpy as np
# import pyxyz.test
# pyxyz.test.run_tests()

if __name__ == "__main__":
    p = Confpool()
    p.include_from_file("crest_conformersA.xyz")
    p.include_from_file("crest_conformersB.xyz")
    print("Number of molecules = {}".format(p.size))

    H2KC = 627.509474063
    KC2H = 1 / H2KC
    n_del = 0
    p["Energy"] = lambda m: float(m.descr.strip())
    n_del += p.upper_cutoff("Energy", 5.0 * KC2H)
    print("{} confs were filtered in total".format(n_del))
    p.sort("Energy") # Ascending by default. Use ascending=False to override
    n_del += p.rmsd_filter(0.3)['DelCount']
    print("{} confs were filtered in total".format(n_del))
    p.descr = lambda m: "Conf #{}; Energy = {:6f} a.u.".format(m.idx + 1, m["Energy"])
    p.save('check.xyz')

    gjftext = """\
%nprocs=24
# opt=(tight,maxcycle=300) freq=noraman def2tzvp empiricaldispersion=gd3bj rpbe1pbe scrf=(cpcm,solvent=water)

Conformation_{idx}

0 1
{xyz}


"""
    for i in range(p.size):
        xyz = p[i].xyz # xyz is a Numpy matrix
        atom_symbols = p.atom_symbols
        xyzlines = []
        for j, symbol in enumerate(atom_symbols):
            xyzlines.append("%2s %14.6f %14.6f %14.6f" % (symbol, *xyz[j]))
        with open("conf_{}.gjf".format(i), 'w') as f:
            f.write(gjftext.format(xyz="\n".join(xyzlines), idx=i))

    p["H2O dist"] = lambda m: m.l(26, 27)
    for i in range(p.size):
        print("Length = {}".format(p[i].l(26, 27)))
        # same as 
        print("Length = {}".format(p[i]["H2O dist"]))
    # Look at the printed numbers. And then
    print("Found {} good structures".format(p.count(lambda m: m["H2O dist"] < 1.8)))
    
    def quality_check(m):
        if m["H2O dist"] < 1.8:
            comment = "Nice"
        else:
            comment = "Garbage"
        return "{}; E={:3f} a.u.".format(comment, m["Energy"])
    for i in range(p.size):
        # Print new descriptions in the console
        print(quality_check(p[i]))
    p.descr = quality_check # Set new descriptions
    p.save('check.xyz')
    
    hbond_condition = lambda m: m.l(25, 10) < 2.0 and m.v(24, 25, 10) > 150.0 and abs(m.z(26, 23, 24, 25)) > 120.0
    num = p.count(hbond_condition)
    print("There are {} matches".format(num))
    p.filter(hbond_condition)

    min_ener = min(p["Energy"]) # p["Energy"] creates a list of floats
    p["E.Rel"] = lambda m: H2KC * (m["Energy"] - min_ener)
    p.descr = lambda m: "{}; E(rel.)={:3f} kcal/mol".format(m.descr, m["E.Rel"])
    p.save('check.xyz')

    import pandas as pd
    p = Confpool()
    p.include_from_file("crest_conformersA.xyz")
    p.include_from_file("crest_conformersB.xyz")
    p["Energy"] = lambda m: float(m.descr.strip())
    p.sort("Energy")

    rmsd_res = p.rmsd_filter(0.3)
    newp = Confpool()
    newp.include_subset(p, [rmsd_res["MinRMSD_pairA"], rmsd_res["MinRMSD_pairB"]])
    newp.save("rmsd_test.xyz")
    print(repr(rmsd_res))

    for i in range(p.size):
        if hbond_condition(p[i]) and p[i].l(26, 27) < 1.8:
            p[i]["MyCondition"] = 1.0
        else:
            p[i]["MyCondition"] = 0.0
    print(pd.DataFrame(p.as_table()))
