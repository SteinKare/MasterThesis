from functions import *
from matplotlib import rcParams
rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
rcParams["font.size"] = "16"
mu = 1
t = 1

labels=[r"$-\pi/a$" ,r"$0$",r"$\pi/a$"]
ticks=[-np.pi, 0, np.pi]

#Mesh
k = np.linspace(-np.pi, np.pi, 100)
Kx, Ky = np.meshgrid(k, k)

#Energy bands
Ep = lp(Kx, Ky, mu, t)
Em = lm(Kx, Ky, mu, t)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_box_aspect(aspect = (1,1,1))
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.set_zlabel("Energy")
surf = ax.plot_surface(Kx, Ky, Ep, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(Kx, Ky, Em, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title("Band Diagram")
plt.show()

#Wave function components
Um1 = um1(Kx, Ky, mu, t)
Um2 = um2(Kx, Ky, mu, t)
Up1 = up1(Kx, Ky, mu, t)
Up2 = up2(Kx, Ky, mu, t)
#Wave functions
Um = np.array([Um1, Um2])
Up = np.array([Up1, Up2])

fig, ax = plt.subplots()
c=ax.imshow((np.angle(Um[1])-np.angle(Um[0])).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
plt.title("Phase difference in lowest-lying eigenvector components")
fig.colorbar(c, ax = ax)
plt.show()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_box_aspect(aspect = (1,1,1))
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
surf = ax.plot_surface(Kx, Ky, np.abs(Um[0]), cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title("Magnitude of first component of lowest-lying eigenvector")
plt.show()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_box_aspect(aspect = (1,1,1))
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
surf = ax.plot_surface(Kx, Ky, np.abs(Um[1]), cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title("Magnitude of second component of lowest-lying eigenvector")
plt.show()

#Derivatives of wave function components
dkx, dky = symbolic_diff_u()

#Berry connections
Am_kx = 1j * (np.conjugate(Um1) * dkx[2](Kx, Ky, mu, t) + np.conjugate(Um2) * dkx[3](Kx, Ky, mu, t))
Am_ky = 1j * (np.conjugate(Um1) * dky[2](Kx, Ky, mu, t) + np.conjugate(Um2) * dky[3](Kx, Ky, mu, t))

fig, ax = plt.subplots()
c=ax.imshow(np.log(np.abs(Am_kx)).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
plt.title("Berry connection Akx of lowest-lying eigenvector")
fig.colorbar(c, ax = ax)
plt.show()

fig, ax = plt.subplots()
c=ax.imshow(np.log(np.abs(Am_ky)).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
plt.title("Berry connection Aky of lowest-lying eigenvector")
fig.colorbar(c, ax = ax)
plt.show()

fig, ax = plt.subplots()
c=ax.imshow(np.angle(Am_kx).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
plt.title("Berry connection Akx of lowest-lying eigenvector")
fig.colorbar(c, ax = ax)
plt.show()

fig, ax = plt.subplots()
c=ax.imshow(np.angle(Am_ky).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
plt.title("Berry connection Aky of lowest-lying eigenvector")
fig.colorbar(c, ax = ax)
plt.show()

#Berry curvatures
Fm, Fp = Berry_curvatures()
Fm_arr = Fm(Kx, Ky, mu, t)
Fp_arr = Fp(Kx, Ky, mu, t)
F1_func, F2_func = Berry_curvature_H_derivative()
F1 = F1_func(Kx, Ky, mu, t)[0, 0]
F2 = F2_func(Kx, Ky, mu, t)[0, 0]

fig, ax = plt.subplots()
c=ax.imshow(np.real(Fm_arr).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.set_xlabel(r"$k_x$")
ax.set_ylabel(r"$k_y$")
clb = plt.colorbar(c)
#clb.ax.set_title(r"$\Omega_{-} (\vec{k})$", y = 1.01)
plt.tight_layout()
plt.savefig(f"BerryCurvature{mu}.pdf", dpi=fig.dpi, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
c=ax.imshow(np.real(F1).T, origin="lower", extent = [np.amin(k), np.amax(k), np.amin(k), np.amax(k)])
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.set_xlabel(r"$k_x$")
ax.set_ylabel(r"$k_y$")
clb = plt.colorbar(c)
#.ax.set_title(r"$F_{-} (\vec{k})$", y = 1.01)
plt.tight_layout()
plt.show()

#Chern numbers
#Check definitions here to see if things are reconcilable

# Cm = 1/(2*np.pi) * np.sum(Fm_arr[~np.isnan(Fm_arr)]) * (k[1] - k[0])**2
# print(Cm)
# Cp = 1/(2*np.pi) * np.sum(Fp_arr[~np.isnan(Fp_arr)]) * (k[1] - k[0])**2
# print(Cp)

# Fm, Fp = Berry_curvatures()
xticks=[-6, -5, -4, -3, -2, -1, 0, 1, 2]
yticks = [-1, 0, 1]
lines = [-4, -2, 0]
mus = np.array([-5.5, -5, -4.5, -3.5, -3, -2.5, -1.5, -1, -0.5, 0.5, 1, 1.5])*t

Cms = np.zeros(mus.shape)
Cps = np.zeros(mus.shape)
for i, mu in enumerate(mus):
    Fm_arr = Fm(Kx, Ky, mu, t)
    Fp_arr = Fp(Kx, Ky, mu, t)
    Cms[i] = 1/(2*np.pi) * np.sum(Fm_arr[~np.isnan(Fm_arr)]) * (k[1] - k[0])**2
    Cps[i] = 1/(2*np.pi) * np.sum(Fp_arr[~np.isnan(Fp_arr)]) * (k[1] - k[0])**2

print(Cms)
print(Cps)
s = 10
c = "black"
fig, ax = plt.subplots()
ax.scatter(mus, Cms, s=s, c=c)
ax.set_xlabel(r"$M$")
ax.set_ylabel(r"$C_{-}$")
for line in lines:
    ax.vlines(line, 0, 1, transform=ax.get_xaxis_transform(), color = "black", linewidth=0.8, linestyle = "solid")
ax.set_xticks(xticks)
ax.set_yticks(yticks)
plt.show()

fig, ax = plt.subplots()
ax.scatter(mus, Cps, s=s, c=c)
ax.set_xlabel(r"$M$")
ax.set_ylabel(r"$C_{+}$")
for line in lines:
    ax.vlines(line, 0, 1, transform=ax.get_xaxis_transform(), color = "black", linewidth=0.8, linestyle = "solid")
ax.set_xticks(xticks)
ax.set_yticks(yticks)
plt.savefig(f"ChernNumberPlus.pdf", dpi=fig.dpi, bbox_inches='tight')
plt.show()

