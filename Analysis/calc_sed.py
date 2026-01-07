import scipy.spatial as sp
import numpy as np
import math


def calc_sed(data):
    coords = data[:,1:4]
    convex_hull = sp.ConvexHull(coords)


    # Area of convex hull; this is a float
    area_hull = convex_hull.area

    # Volume of convex hull; this is a float
    vol_hull = convex_hull.volume

    # Find max distance between vertices for Dmax calculation
    # Necessary for length of corresponding prolate ellipsoid of revolution
    # Note that the "vertices" from ConvexHull are indices, i.e. pointers to
    # the coordinates in "coords" that serve as the vertices.
    vertices = convex_hull.vertices
    Dmax = 0.0
    for r1 in vertices:
        for r2 in vertices:
            dist = math.sqrt((coords[r1][0] - coords[r2][0])**2.0 + \
                             (coords[r1][1] - coords[r2][1])**2.0 + \
                             (coords[r1][2] - coords[r2][2])**2.0 )
            if dist > Dmax:
                Dmax = dist

    #print(area_hull, vol_hull, Dmax)

    pmasses =np.array([
          71.08,
          156.20,
          114.10,
          115.10,
          103.10,
          128.10,
          129.10,
          57.05,
          137.10,
          113.20,
          113.20,
          128.20,
          131.20,
          147.20,
          97.12,
          87.08,
          101.10,
          186.20,
          163.20,
          99.07,
          71.08,
          156.20,
          114.10,
          115.10,
          103.10,
          128.10,
          129.10,
          57.05,
          137.10,
          113.20,
          113.20,
          128.20,
          131.20,
          147.20,
          97.12,
          87.08,
          101.10,
          186.20,
          163.20,
          99.07,
           500,
           62.97])

    vbars=np.array([0.74,
    0.73,
    0.63,
    0.60,
    0.62,
    0.68,
    0.66,
    0.64,
    0.67, 
    0.67,
    0.67,
    0.67,
    0.90,
    0.90,
    0.82,
    0.75,
    0.75,
    0.77,
    0.76,
    0.63,
    0.70,
    0.74,
    0.71,
    0.86,
    0.74,
    0.73,
    0.63,
    0.60,
    0.62,
    0.68,
    0.66,
    0.64,
    0.67, 
    0.67,
    0.67,
    0.67,
    0.90,
    0.90,
    0.82,
    0.75,
    0.75,
    0.77,
    0.76,
    0.63,
    0.70,
    0.74,
    0.71,
    0.86,
    0.65,
    0.501])


    prot_mol_mass = 0
    numerator = 0
    # get total mass:
    for line in data:
        idx = int(line[-1])-1
        prot_mol_mass+=pmasses[idx]
        numerator += (pmasses[idx] * vbars[idx])
    vbar_prot = (numerator / prot_mol_mass) - 0.0025


    # Calculate shell volume from variable of shell thickness
    # Hydration shell thickness is empirically determined to be optimal
    # This expands each hull plane by hydration thickness
    vol_shell_wat = area_hull * 2.8

    # Calculate total volume of hydrated convex hull (including shell waters)
    vol_hyd_hull = vol_hull + vol_shell_wat

    # Estimate axial ratio
    # Minus 3 Ang because DNA duplex a should be rod length, not diagonal
    #   of hull end vertices (Also necessary to make apoferritin axial ratio = 1)
    a = (Dmax/2.0) - 3.0
    b = math.sqrt((3.0 * vol_hull)/(4.0*math.pi*a))

    # But some very spherical structures may not have the diagonal a full 3 Ang longer.
    # Can't have a < b.
    if a > b:
        # Translational Shape Factor
        numerator = math.sqrt(1.0 - (b/a)**2)
        denominator = ((b/a)**0.66666667) * math.log((1 + math.sqrt(1.0 - (b/a)**2))/(b/a))
        Ft = numerator/denominator
    else:
        a = b
        Ft = 1.0

    # Weight shape factor
    # Empirically found to work better with expanded volume
    #  (Many combinations tried)
    Ft = math.sqrt(Ft)

    # Axial ratio of prolate ellipsoid of same volume as convex hull
    a_b_ratio = a/b

    # Find radius of sphere of same volume as hydrated convex hull
    factor = 3.0/(4.0 * math.pi)
    Rht = (factor * vol_hyd_hull)**0.333333
    # Include Shape factor to give effective hydrodynamic translational radius
    Rht = Rht * Ft
    # Rht comes as Angstrom
    # need meters in equation below
    Rh_trans = Rht * 1e-8

    # Calculate Svedberg coeff.
    eta = 0.0100194 # poise
    rho = 0.998234 # g/ml water density
    fT = 6.0*math.pi*eta*Rh_trans
    s = prot_mol_mass * (1.0 - (vbar_prot * rho)) / (6.02214e23 * fT)
    return(s/1e-13)
