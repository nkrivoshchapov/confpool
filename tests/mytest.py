import confpool as m # pymolxyz

# make Confpool < Geompool
if __name__ == "__main__":
    x = m.Confpool()
    x.include("test_res.xyz", energy=lambda x: float(x.strip()))
    # x.energy_filter(1.0, etype='kcal/mol')
    # x.include("testcs.xyz")
    # x.update_description(lambda e, d: "Energy = " + repr(e))
    # x.distance_filter(3, 31, lambda x: x > 1.7)
    x.valence_filter(2, 1, 31, lambda x: x > 160.0)
    x.sort()
    x.save('test_res2.xyz')
