# Use Atlas in place of Create Neighborhoods
def create_neighborhood(r, coords, my_i):
    coord_size = len(coords)
    #print('coordsize ', coord_size)
    #print('r', r)
    list = []
    #ctrlist = []

    #for i in range(coord_size):
    ctr = 0

    for j in range(coord_size):
        if (j != my_i):
            radius = eucl_norm(coords[my_i],coords[j])
            #print('entering first if')
            if (radius<=r):
                #print('entering second if')
                ctr+=1
                list.append(j)

    return list



#Main Loop which iterates over all points and calls different functions
    #for id in range(coord_size):


        #my_list = create_neighborhood(r, coords, id) #Myimplementation
        #my_list = find_neighbors(r,coords) #atlasimplemntation

        #print('mylist', my_list)
        #print(len(my_list))
        #print(coords[my_list[id]])

        #sub_coords = getsubdomain(my_list, coords)
        #print("sub coords:", sub_coords)

        #A = constructa(sub_coords)
        #print('Matrix A:', A)

        cnd = linalg.cond(A)
        #print('Condition number of A ', cnd)

        #InvA = inverta(A)
        #p = projection(coords,my_list)

        #differentiation(sub_coords,InvA,coords[id])


        #print('P Matrix:', p)
        #print('***End of loop ', id, '***')

        #end = time.time()
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