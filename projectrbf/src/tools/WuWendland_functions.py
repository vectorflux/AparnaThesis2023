def wu0(r):
    k = 0
    l = k+2 #floor(d/2) + k + 1
    if (r <= 1):
        phi = (1-r**2)
    else:
        phi = 0
    return phi
def wendland0(r):
    k = 0
    l = k+2 #floor(d/2) + k + 1
    if (r <= 1):
        phi = (1-r)**l
    else:
        phi = 0
    return phi

def wendland1(r):
    k = 1
    l = k+2
    if (r <= 1):
        phi = (((1-r)**(l+1))*((l+1)*r +1))
    else:
        phi = 0
    return phi

def wendland2(r):
    k=2
    l=k+2
    if (r <= 1):
        phi = ((1-r)**(l+2))*(((l**2)+(4*l)+3)*(r**2) +(3*l + 6)*r + 3)
    else:
        phi = 0
    return phi

def wendland3(r):
    k=3
    l=k+2

    if (r <= 1):
        phi = ((1-r)**(l+3))*(((l**3) + 9*(l**2) + 23*l + 15)*(r**3) + (6*(l**2) + 36*l + 45)*(r**2) + (15*l + 45)*r + 15 )
    else:
        phi = 0
    return phi

def wendland0_prime(r): #Derivative of Wendland 0 with k =0 , l = 2

    if (r <= 1):
        phi_p = -2*(1 - r)
    else:
        phi_p = 0
    return phi_p


def wendland1_prime(r): #Derivative of Wendland 1 with k =1 , l = 3

    if (r <= 1):
        #phi_p = -16*(r**3) + 33*(r**2) - 18*r +1
        phi_p = 20*(r**4) - 60*(r**3) + 60*(r**2) - 20*r
    else:
        phi_p = 0
    return phi_p

def wendland_3_1_prime(r): #Derivative of Wendland 1 with k =1 , l = 3

    if (r <= 1):
        phi_p = 20*(r**4) - 60*(r**3) + 60*(r**2) - 20*r 
    else:
        phi_p = 0
    return phi_p

def wendland2_prime(r): #Derivative of Wendland 2 with k =2 , l = 4

    if (r <= 1):
        phi_p = 280*(r**7) - 1344*(r**6) + 2520*(r**5) - 2240*(r**4) + 840*(r**3) - 56*r

    else:
        phi_p = 0
    return phi_p

def wendland3_prime(r): #Derivative of Wendland 3 with k =3 , l = 5

    if (r <= 1):
        phi_p = 5280*(r**10) - 34650*(r**9) + 95040*(r**8) - 138600*(r**7) + 110880*(r**6) - 41580*(r**5) + 3960*(r**3) - 330*r
    else:
        phi_p = 0
    return phi_p
