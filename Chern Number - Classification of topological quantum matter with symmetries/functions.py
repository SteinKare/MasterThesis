import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from matplotlib import rcParams
rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
rcParams["font.size"] = "10"

def R0(k):
    return 0

def R1(k, t):
    return t * np.sin(k)

def R2(k, t):
    return t * np.sin(k)

def R3(kx, ky, mu, t):
    return 2*t + mu - t * (np.cos(kx) + np.cos(ky))

def R(kx, ky, mu, t):
    return np.sqrt(R1(kx, t)**2 + R2(ky, t)**2 + R3(kx, ky, mu, t)**2)

def lp(kx, ky, mu, t):
    return R0(kx) + R(kx, ky, mu, t)

def lm(kx, ky, mu, t):
    return R0(kx) - R(kx, ky, mu, t)

def up1(kx, ky, mu, t):
    front = 1/np.sqrt(2*R(kx, ky, mu, t)*(R(kx, ky, mu, t) - R3(kx, ky, mu, t)))
    return front * (R1(kx, t) - 1j * R2(ky, t))
    
def up2(kx, ky, mu, t):
    front = 1/np.sqrt(2*R(kx, ky, mu, t)*(R(kx, ky, mu, t) - R3(kx, ky, mu, t)))
    return front * (R(kx, ky, mu, t) - R3(kx, ky, mu, t))

def um1(kx, ky, mu, t):
    front = 1/np.sqrt(2*R(kx, ky, mu, t)*(R(kx, ky, mu, t) + R3(kx, ky, mu, t)))
    return front * (R1(kx, t) - 1j * R2(ky, t))

def um2(kx, ky, mu, t):
    front = 1/np.sqrt(2*R(kx, ky, mu, t)*(R(kx, ky, mu, t) + R3(kx, ky, mu, t)))
    return front * (-R(kx, ky, mu, t) - R3(kx, ky, mu, t))

def R0_s(k):
    return 0

def R1_s(k, t):
    return t * sp.sin(k)

def R2_s(k, t):
    return t * sp.sin(k)

def R3_s(kx, ky, mu, t):
    return 2*t + mu - t * (sp.cos(kx) + sp.cos(ky))

def R_s(kx, ky, mu, t):
    return sp.sqrt(R1_s(kx, t)**2 + R2_s(ky, t)**2 + R3_s(kx, ky, mu, t)**2)

def lp_s(kx, ky, mu, t):
    return R0_s(kx) + R_s(kx, ky, mu, t)

def lm_s(kx, ky, mu, t):
    return R0_s(kx) - R_s(kx, ky, mu, t)

def up1_s(kx, ky, mu, t):
    front = 1/sp.sqrt(2*R_s(kx, ky, mu, t)*(R_s(kx, ky, mu, t) - R3_s(kx, ky, mu, t)))
    return front * (R1_s(kx, t) - 1j * R2_s(ky, t))
    
def up2_s(kx, ky, mu, t):
    front = 1/sp.sqrt(2*R_s(kx, ky, mu, t)*(R_s(kx, ky, mu, t) - R3_s(kx, ky, mu, t)))
    return front * (R_s(kx, ky, mu, t) - R3_s(kx, ky, mu, t))

def um1_s(kx, ky, mu, t):
    front = 1/sp.sqrt(2*R_s(kx, ky, mu, t)*(R_s(kx, ky, mu, t) + R3_s(kx, ky, mu, t)))
    return front * (R1_s(kx, t) - 1j * R2_s(ky, t))

def um2_s(kx, ky, mu, t):
    front = 1/sp.sqrt(2*R_s(kx, ky, mu, t)*(R_s(kx, ky, mu, t) + R3_s(kx, ky, mu, t)))
    return front * (-R_s(kx, ky, mu, t) - R3_s(kx, ky, mu, t))
    
def symbolic_diff_u():
    kx, ky, mu, t = sp.symbols("kx, ky, mu, t")
    up1 = up1_s(kx, ky, mu, t)
    up2 = up2_s(kx, ky, mu, t)
    um1 = um1_s(kx, ky, mu, t)
    um2 = um2_s(kx, ky, mu, t)
    dkx_up1_s = sp.lambdify([kx, ky, mu, t], sp.diff(up1, kx))
    dkx_up2_s = sp.lambdify([kx, ky, mu, t], sp.diff(up2, kx))
    dkx_um1_s = sp.lambdify([kx, ky, mu, t], sp.diff(um1, kx))
    dkx_um2_s = sp.lambdify([kx, ky, mu, t], sp.diff(um2, kx))
    
    dky_up1_s = sp.lambdify([kx, ky, mu, t], sp.diff(up1, ky))
    dky_up2_s = sp.lambdify([kx, ky, mu, t], sp.diff(up2, ky))
    dky_um1_s = sp.lambdify([kx, ky, mu, t], sp.diff(um1, ky))
    dky_um2_s = sp.lambdify([kx, ky, mu, t], sp.diff(um2, ky))
    dkx = [dkx_up1_s, dkx_up2_s, dkx_um1_s, dkx_um2_s]
    dky = [dky_up1_s, dky_up2_s, dky_um1_s, dky_um2_s]
    return dkx, dky

def Berry_curvatures():
    kx, ky, mu, t = sp.symbols("kx, ky, mu, t", real = True)
    up1 = up1_s(kx, ky, mu, t)
    up2 = up2_s(kx, ky, mu, t)
    um1 = um1_s(kx, ky, mu, t)
    um2 = um2_s(kx, ky, mu, t)
    
    dkx_up1_s = sp.diff(up1, kx)
    dkx_up2_s = sp.diff(up2, kx)
    dkx_um1_s = sp.diff(um1, kx)
    dkx_um2_s = sp.diff(um2, kx)
    
    dky_up1_s = sp.diff(up1, ky)
    dky_up2_s = sp.diff(up2, ky)
    dky_um1_s = sp.diff(um1, ky)
    dky_um2_s = sp.diff(um2, ky)
    dkx = [dkx_up1_s, dkx_up2_s, dkx_um1_s, dkx_um2_s]
    dky = [dky_up1_s, dky_up2_s, dky_um1_s, dky_um2_s]
    
    Am_kx = sp.I*(sp.conjugate(um1) * dkx_um1_s + sp.conjugate(um2) * dkx_um2_s)
    Am_ky = sp.I*(sp.conjugate(um1) * dky_um1_s + sp.conjugate(um2) * dky_um2_s)
    Ap_kx = sp.I*(sp.conjugate(up1) * dkx_up1_s + sp.conjugate(up2) * dkx_up2_s)
    Ap_ky = sp.I*(sp.conjugate(up1) * dky_up1_s + sp.conjugate(up2) * dky_up2_s)
    
    Fm = sp.diff(Am_ky, kx) - sp.diff(Am_kx, ky)
    Fp = sp.diff(Ap_ky, kx) - sp.diff(Ap_kx, ky)
    Fm = sp.lambdify([kx, ky, mu, t], Fm)
    Fp = sp.lambdify([kx, ky, mu, t], Fp)
    return Fm, Fp

def Berry_curvature_H_derivative():
    kx, ky, mu, t = sp.symbols("kx, ky, mu, t", real = True)
    up1 = up1_s(kx, ky, mu, t)
    up2 = up2_s(kx, ky, mu, t)
    um1 = um1_s(kx, ky, mu, t)
    um2 = um2_s(kx, ky, mu, t)
    T1 = sp.Matrix([um1, um2])
    T2 = sp.Matrix([up1, up2])
    
    d1 = R1_s(kx, t)
    d2 = R2_s(ky, t)
    d3 = R3_s(kx, ky, mu, t)

    Ep = lp_s(kx, ky, mu, t)
    Em = lm_s(kx, ky, mu, t)
    sigma_3 = sp.Matrix([[1, 0], [0, -1]]) #Not used here since fermionic system
    Tk = sp.Matrix([[T1[0], T2[0]], 
                    [T1[1], T2[1]]])
    Tk_dag = sp.transpose(sp.conjugate(Tk))
    H = sp.Matrix([[d3, d1 - sp.I * d2],
                  [d1 + sp.I * d2, d3]])
    dHx = sp.diff(H, kx)
    dHy = sp.diff(H, ky)
    Vx = Tk_dag * dHx * Tk
    Vy = Tk_dag * dHy * Tk
    Omega1 = sp.I * ((sp.transpose(sp.conjugate(T1)) * Vx * T2 * sp.transpose(sp.conjugate(T2)) * Vy * T1) - (sp.transpose(sp.conjugate(T1)) * Vy * T2 * sp.transpose(sp.conjugate(T2)) * Vx * T1)) / (Ep - Em)**2
    Omega2 = sp.I * ((sp.transpose(sp.conjugate(T2)) * Vx * T1 * sp.transpose(sp.conjugate(T1)) * Vy * T2) - (sp.transpose(sp.conjugate(T2)) * Vy * T1 * sp.transpose(sp.conjugate(T1)) * Vx * T2)) / (Ep - Em)**2
    Omega1 = sp.lambdify([kx, ky, mu, t], Omega1)
    Omega2 = sp.lambdify([kx, ky, mu, t], Omega2)
    return Omega1, Omega2
    