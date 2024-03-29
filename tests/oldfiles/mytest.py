import os
import confpool as cp
import pandas as pd
# import pyxyz.test
# pyxyz.test.run_tests()

def func(x):
    # print("x is " + repr(x))
    # print("returns " + repr(x["OH contact"]))
    return x["OH contact"] < 2.0

if __name__ == "__main__":
    curdir = os.path.dirname(os.path.realpath(__file__))
    p = cp.Confpool()
    p.include_from_file(os.path.join(curdir, "crest_conformersA.xyz"))
    p.include_from_file(os.path.join(curdir, "crest_conformersB.xyz"))
    print("Number of molecules = {}".format(p.size))

    H2KC = 627.509474063
    KC2H = 1 / H2KC
    n_del = 0
    p["Energy"] = lambda m: float(m.descr.strip())
    n_del += p.upper_cutoff("Energy", 5.0 * KC2H)
    print("{} confs were filtered in total".format(n_del))
    p.sort("Energy") # Ascending by default. Use ascending=False to override
    n_del += p.rmsd_filter(0.3)
    print("{} confs were filtered in total".format(n_del))
    p.descr = lambda m: "Conf #{}; Energy = {:6f} a.u.".format(m.idx + 1, m["Energy"])
    p.save('check.xyz')

    for i in range(p.size):
        xyz = p[i].xyz # xyz is a Numpy matrix
        atom_symbols = p.atom_symbols
        gjftext = """\
%nprocs=24
# opt=(tight,maxcycle=300) freq=noraman def2tzvp empiricaldispersion=gd3bj rpbe1pbe scrf=(cpcm,solvent=water)

Conformation_{idx}

0 1
{xyz}


"""
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
    p.save('check.xyz') # Now show this file to your supervisor
    
    hbond_condition = lambda m: m.l(25, 10) < 2.0 and m.v(24, 25, 10) > 150.0 and abs(m.z(26, 23, 24, 25)) > 120.0
    num = p.count(hbond_condition)
    print("There are {} matches".format(num))
    p.filter(hbond_condition)

    min_ener = min(p["Energy"]) # p["Energy"] creates a list of floats
    p["E.Rel"] = lambda m: H2KC * (m["Energy"] - min_ener)
    p.descr = lambda m: "{}; E(rel.)={:3f} kcal/mol".format(m.descr, m["E.Rel"])
    p.save('check.xyz')

    import pandas as pd
    p = cp.Confpool()
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
        print(repr(p["MyCondition"]))
    print(pd.DataFrame(p.as_table()))
    raise Exception("pause")


if __name__ == "__main__":
    # Merge several??
    # Confpool<py::object>??
    p = cp.Confpool()

    curdir = os.path.dirname(os.path.realpath(__file__))
    p.include_from_file(os.path.join(curdir, "hydrogen_atom.xyz"))
    p["Energy"] = lambda m: float(m.descr.strip().split('=')[1])
    print("Size = {}".format(p.size))

    p = cp.Confpool()
    p.include_from_file(os.path.join(curdir, "aminoacid_ensemble.xyz"))
    p["Energy"] = lambda m: float(m.descr.strip())
    p["OH contact"] = lambda m: m.l(1, 19)
    p["VA"] = lambda m: m.v(7, 19, 1)
    p["TA"] = lambda m: m.z(8, 6, 7, 19)
    # p.filter(lambda m: abs(m['TA']) > 160.0)
    p.sort("VA") # TODO Ascending/Descending???
    # x.filter(lambda x: x["OH contact"] < 2.0)
    print("A Count = {}".format(p.count(func)))
    print("B Count = {}".format(p.count(lambda mol: mol.l(1, 19) < 2.0)))
    # x.filter(func)
    # x.upper_cutoff("OH contact", 2.0)
    # x.lower_cutoff("OH contact", 2.0)
    # x.update_description(lambda x: "D = {}. Positive energy = {} (old = {})".format(x["OH contact"], -x["Energy"], x.descr().strip()))
    p.descr = lambda m: "TA = {} VA = {}".format(m["TA"], m["VA"])
    p.rmsd_filter(0.3)
    p.save('lolkek.xyz')
    print("XXXSize = {}".format(p.size))

    print(repr(p[0].xyz))
    del p[0]
    del p["VA"]
    print(repr(p[0].xyz))
    print("Table = " + str(pd.DataFrame(p.as_table())))
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
