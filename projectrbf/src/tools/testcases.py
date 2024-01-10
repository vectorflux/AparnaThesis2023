import math
import numpy as np



def test2_fun(lam,th):
    g = 1.5454e-6;
    radius = 1.0
    uu0 = 2 * np.pi * radius / (12 * 86400);
    angle = 0.0 #0.349
    h0 = 2.94e4 / (g*6371000*6371000); #1.9024e-4
    Omega = 7.292e-5;

    coriolis_fun = lambda lam, th: 2 * Omega * (np.sin(th) * np.cos(angle) - np.cos(th) * np.cos(lam) * np.sin(angle))

    h0_fun = lambda lam, th: h0 - (1 / g) * (radius * Omega * uu0 + uu0 ** 2 / 2) * (np.sin(th) * np.cos(angle) - np.cos(lam) * np.cos(th) * np.sin(angle)) ** 2

    lat_fun = uu0 * (np.cos(th) * np.cos(angle) + np.sin(th) * np.cos(lam) * np.sin(angle))
    lon_fun = -uu0 * np.sin(angle) * np.sin(lam)

    h = h0_fun(lam,th)

    f = coriolis_fun(lam,th)

    return lon_fun,lat_fun, h, f


def test6_fun(lam,th):

    g = 1.5454e-6;
    Omega = 7.292e-5;
    omega = 7.848e-6;
    K = omega;
    h0 = 1.25568e-3;#8000/6371000
    R = 4
    radius = 1.0
    uu0 = 2 * np.pi * radius / (12 * 86400);
    angle = 0.0

    coriolis_fun = lambda lam, th: 2 * Omega * (np.sin(th))
    #h0 = 2.94e4 / (g*6371000*6371000); #1.9024e-4

    A_fun = lambda lam, th: 0.5 * omega * (2 * Omega + omega) * np.cos(th) ** 2 + 0.25 * K * K * np.cos(th) ** (2 * R) * ((R + 1) * np.cos(th) ** 2 + (2 * R * R - R - 2) - 2 * R * R * np.cos(th) ** (-2))
    B_fun = lambda lam, th: 2 * (Omega + omega) * K * np.cos(th) ** R * ((R * R + 2 * R + 2) - (R + 1) ** 2 * np.cos(th) ** 2) / ((R + 1) * (R + 2))
    C_fun = lambda lam, th: 0.25 * K * K * np.cos(th) ** (2 * R) * ((R + 1) * np.cos(th) ** 2 - (R + 2))

    h0_fun = lambda lam, th: h0 + radius * (radius / g )* (A_fun(lam, th) + B_fun(lam, th) * np.cos(R * lam) + C_fun(lam, th) * np.cos(2 * R * lam))
    u_fun = lambda lam, th: radius * omega * np.cos(th) + radius * K * np.cos(th) ** (R - 1) * (R * np.sin(th) ** 2 - np.cos(th) ** 2) * (np.cos(R * lam))
    v_fun = lambda lam, th: (-radius * K * R * np.cos(th) ** (R - 1) * np.sin(th) * np.sin(R * lam))


    lat_fun = u_fun(lam,th)
    lon_fun = v_fun(lam,th)
    h = h0_fun(lam,th)
    f = coriolis_fun(lam,th)

   # s2c_lat = np.array([-np.sin(lam), np.cos(lam), 0.])
   # s2c_lon = np.array([-np.cos(lam)*np.sin(th), -np.sin(lam)*np.sin(th), np.cos(th)])

    #uvw = lat_fun*s2c_lat + lon_fun*s2c_lon

    #uvw = cartesian_fromlonlat(lat_fun, lon_fun, lam, th)

    return lon_fun,lat_fun, h, f


def test0_fun(lam,th):

    h0 = 1000/6371000
    h1 = 500/6371000

    a = 1.0 #radius of unit sphere

    R = a/3

    snp = np.sin(np.pi/2)
    csp = np.cos(np.pi/2)

    r = a*np.arccos(snp*np.sin(th) + csp*np.cos(th)*np.cos(lam-(np.pi)))

    if r<R:
        h = (h0/2)*(1 + np.cos((np.pi*r)/(R)))
    else: 
        h = 0.0
    
    lon_fun = 0.0
    lat_fun = 0.0
    f = 0.0

    return lon_fun,lat_fun, h, f
