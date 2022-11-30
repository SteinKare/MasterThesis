import numpy as np
import sympy as sp
from scipy.integrate import dblquad
from scipy.optimize import root

def shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    
       eta_tilde = J1 * D1 / 2 * ( sp.exp(-sp.I * ky) + 2 * sp.exp(sp.I * ky/2) * sp.cos(sp.sqrt(3)/2 * kx) ) 
    
       eta_tilde_squared = J1**2*D1**2/4*(1 + 4*sp.cos(sp.sqrt(3)/2*kx)*sp.cos(3/2*ky) + 4*(sp.cos(sp.sqrt(3)/2*kx))**2) # This value corresponds to |\tilde{\eta}_k|^2 in eq. 5.4.56, and is rewritten in terms of real functions only for the sake of numerical stability
       
       xi = 4*sp.sin(sp.sqrt(3)/2*kx)*(sp.cos(3/2*ky) - sp.cos(sp.sqrt(3)/2*kx))          # xi_k eq. 5.4.58
       
       zeta = 2*(2*sp.cos(sp.sqrt(3)/2*kx)*sp.cos(3/2*ky) + sp.cos(sp.sqrt(3)*kx))          # zeta_k eq. 5.4.59
       
       im_psi = 0.5*J2*D2*xi                                                   # Im(psi_k) in eq. 5.4.56
       
       re_psi = 0.5*G*D3*zeta                                                    # Re(psi_k) in eq. 5.4.56
       
       psi = re_psi + sp.I * im_psi
       
       im_tau = 0.25*D*D2*zeta                                                   # Im(tau_k) in eq. 5.4.56
       
       re_tau = 0.25*D*D3*xi                                                   # Re(tau_k) in eq. 5.4.56
       
       tau = re_tau + sp.I * im_psi
       
       real_psi_tau_conj = sp.re(psi * sp.conjugate(tau))
       
       sqrt_cap_Xi = 2*sp.sqrt(eta_tilde_squared*(re_tau**2 + im_psi**2) + (re_psi*re_tau + im_psi*im_tau)**2)            # This value corresponds to the squareroot of eq. 5.4.79
      
       Ep = sp.sqrt(mu**2 - eta_tilde_squared - im_psi**2 - re_psi**2 - im_tau**2 - re_tau**2 + sqrt_cap_Xi) # plus version of eq. 5.4.80
       
       Em = sp.sqrt(mu**2 - eta_tilde_squared - im_psi**2 - re_psi**2 - im_tau**2 - re_tau**2 - sqrt_cap_Xi) # minus version of eq. 5.4.80
       
       Lambda = tau * sp.conjugate(tau) + psi * sp.conjugate(psi) + (sp.conjugate(tau))**2 - (sp.conjugate(psi))**2
       
       re_tau_minus_i_im_psi = re_tau - sp.I * im_psi
       
       Sigma_p = 2 * eta_tilde_squared * re_tau_minus_i_im_psi + (sp.conjugate(tau) + sp.conjugate(psi)) * (2 * real_psi_tau_conj + sqrt_cap_Xi)
       
       Sigma_m = 2 * eta_tilde_squared * re_tau_minus_i_im_psi + (sp.conjugate(tau) + sp.conjugate(psi)) * (2 * real_psi_tau_conj - sqrt_cap_Xi)
      
       return Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m
   
#Matrix elements

def T11(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m = shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    return (mu + Ep) * (2 * real_psi_tau_conj - sqrt_cap_Xi) / (sp.conjugate(eta_tilde) * (sqrt_cap_Xi - Lambda))

def T21(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m = shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    return 2 * (mu + Ep) * re_tau_minus_i_im_psi / (sqrt_cap_Xi - Lambda)

def T31(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m = shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    return Sigma_m / (sp.conjugate(eta_tilde) * (sqrt_cap_Xi - Lambda))

def T12(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m = shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    return - (mu + Em) * (2 * real_psi_tau_conj + sqrt_cap_Xi) / (sp.conjugate(eta_tilde) * (sqrt_cap_Xi + Lambda))

def T22(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m = shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    return - 2 * (mu + Em) * re_tau_minus_i_im_psi / (sqrt_cap_Xi + Lambda)

def T32(kx, ky, mu, D1, D2, D3, J1, J2, G, D):
    Ep, Em, real_psi_tau_conj, sqrt_cap_Xi, eta_tilde, Lambda, re_tau_minus_i_im_psi, Sigma_p, Sigma_m = shorthands(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    return - Sigma_p / (sp.conjugate(eta_tilde) * (sqrt_cap_Xi + Lambda))

def Berry_curvature_1():
    kx, ky, mu, D1, D2, D3, J1, J2, G, D = sp.symbols("kx, ky, mu, D1, D2, D3, J1, J2, G, D", real = True)
    
    t11 = T11(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    t21 = T21(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    t31 = T31(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    
    dkx_T11 = sp.diff(t11, kx)
    dkx_T21 = sp.diff(t21, kx)
    dkx_T31 = sp.diff(t31, kx)
    
    dkx_T11_conj = sp.diff(sp.conjugate(t11), kx)
    dkx_T21_conj = sp.diff(sp.conjugate(t21), kx)
    dkx_T31_conj = sp.diff(sp.conjugate(t31), kx)
    
    dky_T11 = sp.diff(t11, ky)
    dky_T21 = sp.diff(t21, ky)
    dky_T31 = sp.diff(t31, ky)
    
    dky_T11_conj = sp.diff(sp.conjugate(t11), ky)
    dky_T21_conj = sp.diff(sp.conjugate(t21), ky)
    dky_T31_conj = sp.diff(sp.conjugate(t31), ky)
    
    Omega = sp.I * (dkx_T11_conj * dky_T11 - dky_T11_conj * dkx_T11 + dkx_T21_conj * dky_T21 - dky_T21_conj * dkx_T21 - (dkx_T31_conj * dky_T31 - dky_T31_conj * dkx_T31) )
    return sp.lambdify([kx, ky, mu, D1, D2, D3, J1, J2, G, D], Omega)
    
    
def Berry_curvature_2():
    kx, ky, mu, D1, D2, D3, J1, J2, G, D = sp.symbols("kx, ky, mu, D1, D2, D3, J1, J2, G, D", real = True)
    
    t12 = T12(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    t22 = T22(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    t32 = T32(kx, ky, mu, D1, D2, D3, J1, J2, G, D)
    
    dkx_T12 = sp.diff(t12, kx)
    dkx_T22 = sp.diff(t22, kx)
    dkx_T32 = sp.diff(t32, kx)
    
    dkx_T12_conj = sp.diff(sp.conjugate(t12), kx)
    dkx_T22_conj = sp.diff(sp.conjugate(t22), kx)
    dkx_T32_conj = sp.diff(sp.conjugate(t32), kx)
    
    dky_T12 = sp.diff(t12, ky)
    dky_T22 = sp.diff(t22, ky)
    dky_T32 = sp.diff(t32, ky)
    
    dky_T12_conj = sp.diff(sp.conjugate(t12), ky)
    dky_T22_conj = sp.diff(sp.conjugate(t22), ky)
    dky_T32_conj = sp.diff(sp.conjugate(t32), ky)
    
    Omega = sp.I * (dkx_T12_conj * dky_T12 - dky_T12_conj * dkx_T12 + dkx_T22_conj * dky_T22 - dky_T22_conj * dkx_T22 - (dkx_T32_conj * dky_T32 - dky_T32_conj * dkx_T32) )
    return sp.lambdify([kx, ky, mu, D1, D2, D3, J1, J2, G, D], Omega)
    
def BZ_integration(f, mu, D1, D2, D3, J1, J2, G, D):
    return dblquad(f, -2*np.pi/np.sqrt(3), 2*np.pi/np.sqrt(3), lambda x : 1/np.sqrt(3) * np.abs(x), lambda x : 4*np.pi/3 - 1/np.sqrt(3) * np.abs(x), args=(mu, D1, D2, D3, J1, J2, G, D))[0]

def BZ_kvecs(N):
    #Construct the grid that the sum over the first Brillouin zone is to be taken over
    Kx = np.zeros((N+1, N+1))
    Ky = np.zeros((N+1, N+1))
    for nx in range(0, N+1):
        for ny in range(0, N+1):
            Kx[nx, ny] = (nx - ny)*2*np.pi/(N*np.sqrt(3)) + 1e-6
            Ky[nx, ny] = (nx + ny)*2*np.pi/(3*N) + 1e-6
    return Kx, Ky

def grid_ksum(func, Kx, Ky, D1, D2, D3, mu, J1, J2, G, D, N):
    #Calculate the sum of func over the first Brillouin zone using predefined grid
    #Avoids performing the double for loop more than once
    func_on_grid = func(Kx, Ky, D1, D2, D3, mu, J1, J2, G, D)
    #Remove equivalent BZ points
    pt1 = func(Kx[N,N], Ky[N,N], D1, D2, D3, mu, J1, J2, G, D)
    pt2 = func(Kx[0,N], Ky[0,N], D1, D2, D3, mu, J1, J2, G, D)
    return (np.sum(func_on_grid) - pt1 - pt2)
