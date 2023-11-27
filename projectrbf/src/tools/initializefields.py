import numpy as np
from operator_matrices import *

def initialize_fields(uvwh, n_p, ic_type = "case_6"):


    if ic_type == "case_2":
        g = 9.80616;
        h0 = 2.94e4 / g;
        Omega = 7.292e-5;
        uu0 = 2 * np.pi * radius / (12 * 86400);
        angle = 0
        h0_fun = lambda lam, th: h0 - 1 / g * (radius * Omega * uu0 + uu0 ** 2 / 2) * (
                    np.sin(th) * np.cos(angle) - np.cos(lam) * np.cos(th) * np.sin(angle)) ** 2
        u_fun = lambda lam, th: uu0 * (np.cos(th) * np.cos(angle) + np.sin(th) * np.cos(lam) * np.sin(angle))
        v_fun = lambda lam, th: -uu0 * np.sin(angle) * np.sin(lam)

        coriolis_fun = lambda lam, th: 2 * Omega * (
                    np.sin(th) * np.cos(angle) - np.cos(th) * np.cos(lam) * np.sin(angle))
    elif ic_type == "case_6":
        g = 9.80616;
        Omega = 7.292e-5;
        omega = 7.848e-6;
        K = omega;
        h0 = 8e3;
        R = 4
        A_fun = lambda lam, th: 0.5 * omega * (2 * Omega + omega) * np.cos(th) ** 2 + 0.25 * K * K * np.cos(th) ** (
                    2 * R) * ((R + 1) * np.cos(th) ** 2 + (2 * R * R - R - 2) - 2 * R * R * np.cos(th) ** (-2))
        B_fun = lambda lam, th: 2 * (Omega + omega) * K * np.cos(th) ** R * (
                    (R * R + 2 * R + 2) - (R + 1) ** 2 * np.cos(th) ** 2) / ((R + 1) * (R + 2))
        C_fun = lambda lam, th: 0.25 * K * K * np.cos(th) ** (2 * R) * ((R + 1) * np.cos(th) ** 2 - (R + 2))

        h0_fun = lambda lam, th: h0 + radius * radius / g * (
                    A_fun(lam, th) + B_fun(lam, th) * np.cos(R * lam) + C_fun(lam, th) * np.cos(2 * R * lam))
        u_fun = lambda lam, th: radius * omega * np.cos(th) + radius * K * np.cos(th) ** (R - 1) * (
                    R * np.sin(th) ** 2 - np.cos(th) ** 2) * (np.cos(R * lam))
        v_fun = lambda lam, th: -radius * K * R * np.cos(th) ** (R - 1) * np.sin(th) * np.sin(R * lam)

        coriolis_fun = lambda lam, th: 2 * Omega * np.sin(th)
    else:
        raise Exception(f"Unknown initial condition type: {ic_type}")


    for i in range(n_p):
            local_pos_x = x_c[i] + 0.5 * hx * unif2d_x
            local_pos_y = y_c[j] + 0.5 * hy * unif2d_y
            local_qp_x = x_c[i] + 0.5 * hx * pts2d_x
            local_qp_y = y_c[j] + 0.5 * hy * pts2d_y
            h[i, j, 0, :] = h0_fun(local_pos_x, local_pos_y)
            u[i, j, 0, :] = u_fun(local_pos_x, local_pos_y)
            v[i, j, 0, :] = v_fun(local_pos_x, local_pos_y)
            coriolis[i, j, 0, :] = coriolis_fun(local_qp_x, local_qp_y)


#def construct_rhs(uvwh, allD, allP, xyz ,nrj_size_list, n_p):


 #   uvwh0, f = initialize_fields()

#  RHS_D = construct_rhsd()

 #   Ru = -px*(RHS_D)

  #  Rv = -py*(RHS_D)

   # Rw = -pz*(RHS_D)

    #Rh =

    #return Ru, Rv, Rw, Rh


#def construct_rhsd():
def construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD):
    n_p0 = len(nrj_size_list)
    rhsd = np.empty(3)
    it = 0
    tmp_xyz = np.empty(3)
    termA = np.empty(3)
    termB = np.empty(3)
    termC = np.empty(3)
    termCx = np.empty(n_p0)
    termCy = np.empty(n_p0)
    termCz = np.empty(n_p0)

    Ru = np.empty(n_p0)
    Rv = np.empty(n_p0)
    Rw = np.empty(n_p0)
    Rh = np.empty(n_p0)

    g = 9.80616
    Omega = 7.292e-5
    f = 1.4e-3

    #need to get neighborhood Dnx, Dny, Dnz (saved sequentially for each point)


    def get_Dnxyz_from_allD(allD, it, k):

         #gives the length of nrj for the index point id
        
        itx = int(it+k)
        
        Dnx = allD[it: itx, 0]
        Dny = allD[it: itx, 1]
        Dnz = allD[it: itx, 2]
        nearest = allnearest[it:itx]

        return Dnx, Dny, Dnz, nearest

    # need to get uvwh (saved at each point)
    def get_uvwh_r(uvwh,nearest):
        uvwh_r = np.empty([len(nearest),4])
        for n in range(len(nearest)):
            uvwh_r[n] = uvwh[int(nearest[n]),:]

        return uvwh_r

    def getpxyz(xyz):

        px = np.ndarray(3)
        py = np.ndarray(3)
        pz = np.ndarray(3)

        x = xyz[0]
        y = xyz[1]
        z = xyz[2]

        px[0] = (1 - (x ** 2))
        px[1] = (0 - (x * y))
        px[2] = (0 - (x * z))

        py[0] = (0 - (x * y))
        py[1] = (1 - (y ** 2))
        py[2] = (0 - (y * z))

        pz[0] = (0 - (x * z))
        pz[1] = (0 - (y * z))
        pz[2] = (1 - (z ** 2))


        return px[None], py[None], pz[None]


    for i in range(n_p0):

        k = int(nrj_size_list[i])
        Dnx, Dny, Dnz, nearest = get_Dnxyz_from_allD(allD, it,k)
        uvwh_r = get_uvwh_r(uvwh, nearest)

        u = uvwh_r[:,0]
        v = uvwh_r[:, 1]
        w = uvwh_r[:, 2]
        h = uvwh_r[:, 3]

        # u.(Dnx.u_r) + v.(Dny.u_r) + w.(Dnz. u_r)
        termA[0] = (uvwh[i,0]*(np.dot(np.transpose(Dnx),u)) +
                 uvwh[i,1]*(np.dot(np.transpose(Dny),u)) +
                 uvwh[i,2]*(np.dot(np.transpose(Dnz),u)))

        termA[1] = (uvwh[i, 0] * (np.dot(np.transpose(Dnx), v)) +
                       uvwh[i, 1] * (np.dot(np.transpose(Dny), v)) +
                       uvwh[i, 2] * (np.dot(np.transpose(Dnz), v)))

        termA[2] = (uvwh[i, 0] * (np.dot(np.transpose(Dnx), w)) +
                    uvwh[i, 1] * (np.dot(np.transpose(Dny), w)) +
                    uvwh[i, 2] * (np.dot(np.transpose(Dnz), w)))

        it = it + k

        #f. [3,1]
        termB[0] = xyz[i,1]*uvwh[i,2] - xyz[i,2] *uvwh[i,1]
        termB[1] = xyz[i,2]*uvwh[i,0] - xyz[i,0] *uvwh[i,2]
        termB[2] = xyz[i,0]*uvwh[i,1] - xyz[i,1] *uvwh[i,0]

        #g.[Dnx Dny Dnz]. h

        termCx= Dnx
        termCy= Dny
        termCz= Dnz
        
        term_c = np.row_stack([termCx,termCy,termCz])
        termC = np.dot(term_c,h)

        rhsd[0] = termA[0] + f*(termB[0]) + g*(termC[0])
        rhsd[1] = termA[1] + f*(termB[1]) + g*(termC[1])
        rhsd[2] = termA[2] + f*(termB[2]) + g*(termC[2])

        rhsd[:,None]


        #Get px, py, pz at the node where we are calculating
        px, py, pz = getpxyz(xyz[i])

        Rh[i] = ( uvwh[i, 0]*(np.dot(Dnx,h))   + uvwh[i, 1]*(np.dot(Dny,h)) + uvwh[i, 2]*(np.dot(Dnz,h)) +
            uvwh[i,3]*(np.dot(Dnx,u) + np.dot(Dny,v)+ np.dot(Dnz,w)))

        Ru[i] = -np.dot(px, rhsd)
        Rv[i] = -np.dot(py, rhsd)
        Rw[i] = -np.dot(pz, rhsd)



    return Rh, Ru, Rv, Rw

def set_initial_conditions(uvwh, xyz, n_p, ghost):

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    for n in range(n_p):
        if not ghost[n]: #initializing all internal values
            uvwh[n,0] = (x[n]**2 +y[n]**2 +z[n]**2)
            uvwh[n,1] = 2*(x[n]**2 +y[n]**2 +z[n]**2)
            uvwh[n,2] = 3*(x[n]**2 +y[n]**2 +z[n]**2)
            uvwh[n,3] = 4*(x[n]**2 +y[n]**2 +z[n]**2)
        else:
            uvwh[n,0] = 0
            uvwh[n,1] = 0
            uvwh[n,2] = 0
            uvwh[n,3] = 0


    return uvwh

def validate_halo_exchange(uvwh,xyz, n_p,ghost):
    vt = 1
    for n in range(n_p):
        if ghost[n]:
            if (uvwh[n,3] != 4*(xyz[n,0]**2 +xyz[n,1]**2 +xyz[n,2]**2)):
                print("Halo exchange failed")
                vt = 0

    return vt
