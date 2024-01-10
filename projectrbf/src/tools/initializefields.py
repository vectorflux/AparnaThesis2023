import numpy as np
from operator_matrices import *
import math
from helperfuncs import *
from testcases import *


def construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD, ghost):
    n_p = len(uvwh)
    rhsd = np.zeros(3)
    it = 0
    termA = np.zeros(3)
    termB = np.zeros(3)
    termC = np.zeros(3)

    termc_collect = np.zeros((n_p,3))
    #termCx = np.zeros(n_p, dtype = np.float64)
    #termCy = np.zeros(n_p, dtype = np.float64)
    #termCz = np.zeros(n_p, dtype = np.float64)

    Ru = np.zeros(n_p)
    Rv = np.zeros(n_p)
    Rw = np.zeros(n_p)
    Rh = np.zeros(n_p)

    #g = 9.80616
    #f = 1.4e-3
    g = 1.5454e-6

    #need to get neighborhood Dnx, Dny, Dnz (saved sequentially for each point)
    def get_Dnxyz_from_allD(allD, it, k):

         #gives the length of nrj for the index point id
        
        itx = it+k
        
        Dnx = allD[it: itx, 0]
        Dny = allD[it: itx, 1]
        Dnz = allD[it: itx, 2]
        nearest = allnearest[it:itx]
        
        
        #print(len(Dnx),k)
        #if it ==0:
        #print("Dnx:",Dnx)
        #print("Dny:",Dny)
        #print("Dnz:",Dnz)
        
        return Dnx[:,np.newaxis], Dny[:,np.newaxis], Dnz[:,np.newaxis], nearest

    # need to get neighborhood uvwh (saved at each point)
    def get_uvwh_r(uvwh,nearest):
        uvwh_r = np.zeros([len(nearest),4])
        for n in range(len(nearest)):
            uvwh_r[n,0] = uvwh[int(nearest[n]),0]
            uvwh_r[n,1] = uvwh[int(nearest[n]),1]
            uvwh_r[n,2] = uvwh[int(nearest[n]),2]
            uvwh_r[n,3] = uvwh[int(nearest[n]),3]
        
        #print(np.where((np.sum(uvwh_r, axis=1))==0)[0])


        return uvwh_r

    def getpxyz(xyz, i):

        px = np.zeros(3)
        py = np.zeros(3)
        pz = np.zeros(3)

        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        

        #if i==1:
            #print("xyz", x,y,z)

        px[0] = (1 - (x ** 2))
        px[1] = (0 - (x * y))
        px[2] = (0 - (x * z))

        py[0] = (0 - (x * y))
        py[1] = (1 - (y ** 2))
        py[2] = (0 - (y * z))

        pz[0] = (0 - (x * z))
        pz[1] = (0 - (y * z))
        pz[2] = (1 - (z ** 2))

        px = px[np.newaxis,:] #row matrix
        py = py[np.newaxis,:]
        pz = pz[np.newaxis,:]

        return px, py, pz
    
    l = 0
    f = uvwh[:,4]
    #for i in range(len(f)):
        
    #print("Coriolis function range:", max(f), min(f))
    #print("U range:", max(uvwh[:,0]), min(uvwh[:,0]))
    #print("v range:", max(uvwh[:,1]), min(uvwh[:,1]))
    #print("w range:", max(uvwh[:,2]), min(uvwh[:,2]))
    #print("h range:", max(uvwh[:,3]), min(uvwh[:,3]))
    #result = np.where(~uvwh.any(axis=1))[1]
    #result = np.where(np.all(np.isclose(uvwh[:], 0), axis=1))
    #print("Rows with non zero:", result)
    #print(np.where((np.sum(uvwh, axis=1))==0.0)[0])
    
    # Calculating RHS values for each of the internal points
    for i in range(n_p):
        if not ghost[i]:

            k = int(nrj_size_list[l]) #Check the requirement of int typecasting here
            l= l+1
            Dnx, Dny, Dnz, nearest = get_Dnxyz_from_allD(allD, it,k)
            
            it = it + k
            
            xyz_r = getneighcoords(nearest, xyz)
            uvwh_r = get_uvwh_r(uvwh, nearest)

            u = uvwh_r[:, 0]
            v = uvwh_r[:, 1]
            w = uvwh_r[:, 2]
            h = uvwh_r[:, 3]

            #for m in range(k):
                #diffh = h[m] - uvwh[int(nearest[m]),3] 
                #dist = eucl_norm(xyz[i,:], xyz_r[m,:])

            #print("Dot product of gradient operator and u",np.dot(np.transpose(Dnx),u))
            #dnxu = np.zeros(3)

            #dnxu[0] = g*np.dot(np.transpose(Dnx),h)
            #dnxu[1] = g*np.dot(np.transpose(Dny),h)
            #dnxu[2] = g*np.dot(np.transpose(Dnz),h)

            #for i in range(len(xyz_r)):
                #product = np.dot(xyz_r[i,:], dnxu)
                #if not np.isclose(product,0.0):
                    #print( " dotproduct:", product)
            #print("Dnxyz:", np.max(dnxu), np.min(dnxu))

            termA[0] = (uvwh[i,0]*(np.dot(np.transpose(Dnx),u)) +
                 uvwh[i,1]*(np.dot(np.transpose(Dny),u)) +
                 uvwh[i,2]*(np.dot(np.transpose(Dnz),u)))

            termA[1] = (uvwh[i, 0] * (np.dot(np.transpose(Dnx), v)) +
                       uvwh[i, 1] * (np.dot(np.transpose(Dny), v)) +
                       uvwh[i, 2] * (np.dot(np.transpose(Dnz), v)))

            termA[2] = (uvwh[i, 0] * (np.dot(np.transpose(Dnx), w)) +
                    uvwh[i, 1] * (np.dot(np.transpose(Dny), w)) +
                    uvwh[i, 2] * (np.dot(np.transpose(Dnz), w)))

            #it = it + k

            #print("Shape of Term A:", np.shape(termA))

            #f. [3,1]
            termB[0] = f[i]*((xyz[i,1]*uvwh[i,2]) - (xyz[i,2]*uvwh[i,1]))
            termB[1] = f[i]*((xyz[i,2]*uvwh[i,0]) - (xyz[i,0]*uvwh[i,2]))
            termB[2] = f[i]*((xyz[i,0]*uvwh[i,1]) - (xyz[i,1]*uvwh[i,0]))

            #g.[Dnx Dny Dnz]. h

            #print("Shape of Term B:", np.shape(termB))

            termCx= np.transpose(Dnx)
            termCy= np.transpose(Dny)
            termCz= np.transpose(Dnz)
            
            #ana_c = 

            term_c = np.row_stack([termCx,termCy,termCz])
            termC = g*(np.dot(term_c,h))

            termc_collect[i,:] = np.transpose(termC)

            #print("shape of Term C:", np.shape(termC)
            #print("TermC:", np.max(termC), np.min(termC), termC)

            #if i==0:
                #print("Shape of termC:", np.shape(termC))
                #print("Shape of termB:", np.shape(termB))
                #print("Shape of termA:", np.shape(termA))

            rhsd[0] = termA[0] + (termB[0]) + (termC[0])
            rhsd[1] = termA[1] + (termB[1]) + (termC[1])
            rhsd[2] = termA[2] + (termB[2]) + (termC[2])

            #rhsd[:,np.newaxis]

            
            #Get px, py, pz at the node where we are calculating
            px, py, pz = getpxyz(xyz[i,:], i)


            Rh[i] = uvwh[i, 0]*(np.dot(np.transpose(Dnx),h)) + uvwh[i, 1]*(np.dot(np.transpose(Dny),h)) + uvwh[i, 2]*(np.dot(np.transpose(Dnz),h)) + uvwh[i, 3]*(np.dot(np.transpose(Dnx),u) + np.dot(np.transpose(Dny),v)+ np.dot(np.transpose(Dnz),w))

            Ru[i] = -np.dot(px, rhsd)
            Rv[i] = -np.dot(py, rhsd)
            Rw[i] = -np.dot(pz, rhsd)
            

            ###################
            

    #print("TermC:", np.max(termc_collect[:,:]), np.min(termc_collect[:,:]))

    return Ru, Rv, Rw, Rh

def set_initial_conditions(xyz, n_p, ghost,lonlat):
    
    N = n_p #internal+halo , 5 field columns
    f = np.zeros(N)

    fields = np.zeros((N,5))
    #Omega = 7.292e-5;
    #angle = 0.0
       
    #lam = np.zeros(len(xyz))
    #th = np.zeros(len(xyz))

    lam = np.radians(lonlat[:,0])
    th = np.radians(lonlat[:,1])

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    
    v_lon = np.zeros(N)
    v_lat = np.zeros(N)
    h = np.zeros(N)

    for n in range(N):
        if not ghost[n]: #calculating all internal values     
            v_lon[n], v_lat[n], h[n], f[n] = test0_fun(lam[n],th[n]) 

        else:
            v_lon[n] = v_lat[n] = h[n] = 0.0
            f[n] = 0.0

            #print("vlon and vlat: ", v_lon[n,0], v_lat[n,0])

    uvw = gvec2cvec_sphere(v_lon, v_lat, lam, th,ghost)
            #uvw, h = test6_fun(lam[n],th[n])
            
    for n in range(n_p):
        if not ghost[n]:
            fields[n,0] = uvw[n,0]
            fields[n,1] = uvw[n,1]
            fields[n,2] = uvw[n,2]
            fields[n,3] = h[n]
            fields[n,4] = f[n]

        else:
            fields[n,0] = 0.
            fields[n,1] = 0.
            fields[n,2] = 0.
            fields[n,3] = 0.
            fields[n,4] = 0.

        
    #print("uvwh and f: ", uvwh)

    return fields


