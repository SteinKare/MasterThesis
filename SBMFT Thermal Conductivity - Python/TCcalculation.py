from functions import *

B1, B2 = Berry_curvature_1(), Berry_curvature_2()

N = 1000
Kx, Ky = BZ_kvecs(N)
dkx = np.abs(Kx[0, 0]-Kx[0, 1])
dky = np.abs(Ky[0, 0]-Ky[0, 1])

J1 = 1
J2s = [5, 2, 0.2]
Gs = [2, 5, 0.2]
D = 0

#Load sols from file:
sols_file = input('Specify file name for solutions matrix: ')
with open(sols_file, "rb") as f:
    sols = np.load(f)

print(sols)
    

for i in range(len(sols)):
    mu = sols[i, 0]
    D1 = sols[i, 1]
    D2 = sols[i, 2]
    D3 = sols[i, 3]
    #B1_int = 1/(2*np.pi)**2 * BZ_integration(B1, mu, D1, D2, D3, J1, J2, G, D)
    #B2_int = 1/(2*np.pi)**2 * BZ_integration(B2, mu, D1, D2, D3, J1, J2, G, D)
    B1_int = 1/(2*np.pi) * grid_ksum(B1, Kx, Ky, D1, D2, D3, mu, J1, J2, G, D, N) * dkx * dky
    B2_int = 1/(2*np.pi) * grid_ksum(B2, Kx, Ky, D1, D2, D3, mu, J1, J2, G, D, N) * dkx * dky
    print("B1:", B1_int)
    print("B2:", B2_int)
    print("Sum", B1_int + B2_int)



    
