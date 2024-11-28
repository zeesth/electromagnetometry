import matplotlib.pyplot as plt
import numpy as np
from dircos import dircos

CM = 1e-7
TnT = 1e9

def distance(size, prism_xyz1, prism_xyz2, mag_incl, mag_decl, field_incl, field_decl, x_theta):
    obs_x = np.arange(start=0, stop=(size+1), step=10, dtype=float)
    obs_y = np.arange(start=0, stop=(size+1), step=10, dtype=float)
    obs_x, obs_y = np.meshgrid(obs_x, obs_y)
    obs_z = np.zeros_like(obs_x)

    obs = [obs_x, obs_y, obs_z]

    prism_x1, prism_y1, prism_z1 = prism_xyz1
    prism_x2, prism_y2 = prism_xyz2[:2]
    
    mag_a, mag_b, mag_c = dircos(mag_incl, mag_decl, x_theta)
    field_a, field_b, field_c = dircos(field_incl, field_decl, x_theta)

    fm1 = mag_a * field_b + mag_b * field_a
    fm2 = mag_a * field_c + mag_c * field_a
    fm3 = mag_a * field_c + mag_c * field_b
    fm4 = mag_a * field_a
    fm5 = mag_b * field_b
    fm6 = mag_c * field_c

    fm = [fm1, fm2, fm3, fm4, fm5, fm6]

    alpha1 = prism_x1 - obs[0]
    alpha2 = prism_x2 - obs[0]
    alpha = [alpha1, alpha2]

    beta1 = prism_y1 - obs[1]
    beta2 = prism_y2 - obs[1]
    beta = [beta1, beta2]

    h = prism_z1 - obs[2]

    return fm, beta, alpha, h, obs

def total_field(mag_it, fm, beta, alpha, h):
    t = 0
    for i in range(2):
        alpha_t = alpha[i]
        alpha_sq = alpha_t ** 2
        for j in range(1, 2):
            beta_t = beta[j]
            beta_sq = beta_t ** 2
            sign = 1

            if i != j:
                sign = -1
            
            r0_sq = alpha_sq + beta_sq + h ** 2
            r0 = np.sqrt(r0_sq)
            r0h = r0 * h

            alphabeta = alpha_t * beta_t

            arg1 = (r0 - alpha_t)/(r0 + alpha_t)
            arg2 = (r0 - beta_t)/(r0 + beta_t)
            arg3 = alpha_sq + r0h + h ** 2
            arg4 = r0_sq + r0h - alpha_sq

            tlog = (fm[2] * np.log(arg1) / 2) + (fm[3] * np.log(arg2) / 2) - (fm[0] * np.log(r0 + h))

            tatan = -(fm[3] * np.arctan2(alphabeta, arg3)) -(fm[4] * np.arctan2(alphabeta, arg4)) +(fm[5] * np.arctan2(alphabeta, r0h))

    t = t + sign * (tlog + tatan)
    t = t * mag_it * CM * TnT

    return t

def user_input(first=True):
    x1 = 1000 * float(input("Prism's x1 (km): "))
    y1 = 1000 * float(input("Prism's y1 (km): "))
    z1 = 1000 * float(input("Prism's z1 (km): "))
    x2 = 1000 * float(input("Prism's x2 (km): "))
    y2 = 1000 * float(input("Prism's y2 (km): "))
    z2 = 1000 * float(input("Prism's z2 (km): "))

    prism_xyz1 = [x1, y1, z1]
    prism_xyz2 = [x2, y2, z2]

    mag_ic = float(input("Magnetization inclination(deg): "))
    mag_d = float(input("Magnetization declination(deg): "))
    mag_it = float(input("Magnetization intensity(A/m): "))

    field_i = float(input("Ambient field inclination(deg): "))
    field_d = float(input("Ambient field declination(deg): "))

    x_theta = float(input("Declination on x-axis(deg): "))
    if first:
        size = 1000 * float(input("Grid size (km): "))
        return size, prism_xyz1, prism_xyz2, mag_ic, mag_d, mag_it, field_i, field_d, x_theta    

    return prism_xyz1, prism_xyz2, mag_ic, mag_d, mag_it, field_i, field_d, x_theta

def plots(obs_x, obs_y, t, t1, t2):
    fig = plt.figure(figsize=[12.8, 4.8])
    y_filter = np.where(obs_y==0)[0][0]

    trid = fig.add_subplot(1, 2, 1, projection="3d",
    xlabel="x (km)",
    ylabel="y (km)",
    zlabel="T (nT)",
    title="Total field anomaly along x")
    trid.plot_surface(obs_x/1000, obs_y/1000, t, cmap="viridis")
    trid.view_init(azim=230)

    twod = fig.add_subplot(1, 2, 2,
    xlabel="x (km)",
    ylabel="T (nT)",
    title="Total field anomaly along x")
    twod.plot(obs_x[0]/1000, t[y_filter, :], label="Sum")
    twod.plot(obs_x[0]/1000, t1[y_filter, :], label="Prism 1")
    twod.plot(obs_x[0]/1000, t2[y_filter, :], label="Prism 2")
    twod.grid()
    twod.set_xlim(0, obs_x[0][-1]/1000)
    twod.legend()
    plt.show()

def main():
    size, prism_xyz1, prism_xyz2, mag_ic, mag_d, mag_it1, field_i1, field_d1, x_theta1 = user_input()

    prism_xyz3, prism_xyz4, mag_ic2, mag_d2, mag_it2, field_i2, field_d2, x_theta2 = user_input(False)

    fm1, beta1, alpha1, h1, obs = distance(size, prism_xyz1, prism_xyz2, mag_ic, mag_d, field_i1, field_d1, x_theta1)
    fm2, beta2, alpha2, h2 = distance(size, prism_xyz2, prism_xyz1, mag_ic, mag_d, field_i1, field_d1, x_theta1)[:-1]

    fm3, beta3, alpha3, h3 = distance(size, prism_xyz3, prism_xyz4, mag_ic2, mag_d2, field_i2, field_d2, x_theta2)[:-1]
    fm4, beta4, alpha4, h4 = distance(size, prism_xyz4, prism_xyz3, mag_ic2, mag_d2, field_i2, field_d2, x_theta2)[:-1]

    t1 = total_field(mag_it1, fm1, beta1, alpha1, h1)
    t2 = total_field(mag_it1, fm2, beta2, alpha2, h2)

    t3 = total_field(mag_it2, fm3, beta3, alpha3, h3)
    t4 = total_field(mag_it2, fm4, beta4, alpha4, h4)

    ta = t1 - t2
    tb = t3 - t4
    tt = ta + tb

    plots(obs[0], obs[1], tt, ta, tb)

main()