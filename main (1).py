import math
import numpy as np
import matplotlib.pyplot as plt
MU = 398600.4415
MUM = 4902.8000
W0 = 7.29 * math.pow(10, -5)
MUS = 132712440018
JD = 0.0
S = 0.0
x = [i * 0 for i in range(6)]
global x1
x1 = [i * 0 for i in range(6)]
global r_r
global m
lum = 0.0
#################################################################################
def JulyData(t,JD):
    year = 2022
    month = 5
    day = 17
    hours = 7
    min = 9
    sec = 33 + t
    a = (14 - month) / 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    JDN = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
    JD = JDN + (hours - 12.) / 24. + min / 1440. + (sec + t) / 86400.0
    return JD

#################################################################################
def TimeStar(t, dd):
    hours = 7
    min = 9
    sec = 33 + t
    v = 0.0027379093
    T = (dd - 2415020.0) / 36525.0
    S0 = 101.25228375 * math.pi / 180 + 360000.77005361 * math.pi / 180 * T + 3.879333333333e-4 * math.pi / 180 * T * T
    M = (sec / 60 + min) / 60 + hours
    S = S0 + M * (1 + v)
    return S
#################################################################################


def GESK(x, Xg, lum):
    Xg[0] = math.cos(S) * x[0] - math.sin(S) * x[1]
    Xg[1] = math.cos(S) * x[0] + math.sin(S) * x[1]
    Xg[2] = x[2]
    lum = math.acos(Xg[0] / math.sqrt(Xg[0] * Xg[0] + Xg[1] * Xg[1]))
    return x, Xg, lum
#################################################################################

def C(x,Ss,xc20,xc30,xc40,xc50,xc60,xc22):
    ae = 6378.136
    J2 = 1082.636023e-6
    J3 = -2.532435e-6
    J4 = -2.37089e-6
    J5 = -0.227716e-6
    J6 = -0.00608e-6
    r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    C22 = 241.29e+8;
    D22 = -136.41e+8;
    J22 = math.sqrt(C22 * C22 + D22 * D22)
    lam22 = math.atan(D22 / C22)
    Re = 6378.155
    A22 = 3 * J22 * MU * Re * Re / math.pow(r, 5)
    lam = math.acos((x[0] * math.cos(Ss) + x[1] * math.sin(Ss)) / math.sqrt(x[0] * x[0] + x[1] * x[1]))

    xc20[0] = 3 / 2 * J2 * MU * ae * ae / math.pow(r, 5) * (5 * x[0] * x[2] * x[2] / r / r - x[0])
    xc20[1] = 3 / 2 * J2 * MU * ae * ae / math.pow(r, 5) * (5 * x[1] * x[2] * x[2] / r / r - x[1])
    xc20[2] = 3 / 2 * J2 * MU * ae * ae / math.pow(r, 5) * (5 * x[2] * x[2] * x[2] / r / r - 3 * x[2])

    xc30[0] = 5 / 2 * J3 * MU * ae * ae * ae / math.pow(r, 7) * (7 * x[0] * x[2] * x[2] * x[2] / r / r - 3 * x[0] * x[2]);
    xc30[1] = 5 / 2 * J3 * MU * ae * ae * ae / math.pow(r, 7) * (7 * x[1] * x[2] * x[2] * x[2] / r / r - 3 * x[1] * x[2]);
    xc30[2] = 5 / 2 * J3 * MU * ae * ae * ae / math.pow(r, 7) * (7 * x[2] * x[2] * x[2] * x[2] / r / r - 3 * x[2] * x[2] + 3 / 5 * r * r);

    xc40[0] = 1 / 4 * J4 * MU * ae * ae * ae * ae / math.pow(r, 7) * x[0] * (315 * math.pow(x[2], 4) / (2 * math.pow(r, 4)) - 105 * math.pow(x[2], 2) /math.pow(r, 2) + 15.0 / 2)
    xc40[1] = 1 / 4 * J4 * MU * ae * ae * ae * ae / math.pow(r, 7) * x[1] * (315 * math.pow(x[2], 4) / (2 * math.pow(r, 4)) - 105 * math.pow(x[2], 2) / math.pow(r, 2) + 15.0 / 2)
    xc40[2] = 1 / 4 * J4 * MU * ae * ae * ae * ae / math.pow(r, 7) * x[2] * (315 * math.pow(x[2], 4) / (2 * math.pow(r, 4)) - 175 * math.pow(x[2], 2) / math.pow(r, 2) + 75.0 / 2)

    xc50[0] = 1 / 8 * J5 * MU * ae * ae * ae * ae * ae / math.pow(r, 9) * x[0] * (693 * math.pow(x[2], 5) / (2 * math.pow(r, 4)) - 630 * math.pow(x[2], 3) / math.pow(r, 2) + 105 * x[2])
    xc50[1] = 1 / 8 * J5 * MU * ae * ae * ae * ae * ae / math.pow(r, 9) * x[1] * (693 * math.pow(x[2], 5) / (2 * math.pow(r, 4)) - 630 * math.pow(x[2], 3) / math.pow(r, 2) + 105 * x[2])
    xc50[2] = 1 / 8 * J5 * MU * ae * ae * ae * ae * ae / math.pow(r, 9) * x[2] * ((693 * math.pow(x[2], 5) - 945 * math.pow(x[2], 3) * r * r) / (2 * math.pow(r, 4)) - 315 * x[2] / math.pow(r,2) - 15 * r * r /x[2])

    xc60[0] = 1 / 16 * J6 * MU * ae * ae * ae * ae * ae * ae / math.pow(r, 9) * x[0] * (3003 * math.pow(x[2], 6) / math.pow(r, 6) - 3465 * math.pow(x[2], 4) / (math.pow(r, 4)) + 945 * math.pow(x[2], 2) / math.pow(r, 2) - 35)
    xc60[1] = 1 / 16 * J6 * MU * ae * ae * ae * ae * ae * ae / math.pow(r, 9) * x[1] * (3003 * math.pow(x[2], 6) / math.pow(r, 6) - 3465 * math.pow(x[2], 4) / (math.pow(r, 4)) + 945 * math.pow(x[2], 2) / math.pow(r,2) - 35)
    xc60[2] = 1 / 16 * J6 * MU * ae * ae * ae * ae * ae * ae / math.pow(r, 9) * x[2] * (3003 * math.pow(x[2], 6) / math.pow(r, 6) - 4851 * math.pow(x[2], 4) / (math.pow(r, 4)) + 2205 * math.pow(x[2], 2) / math.pow(r,2) - 245)

    xc22[0] = A22 * (x[0] * (5 * x[2] * x[2] / r / r - 3) * math.cos(2 * (lum - lam22)) + 2 * x[1] * math.sin(2 * (lum - lam22)))
    xc22[1] = A22 * (x[1] * (5 * x[2] * x[2] / r / r - 3) * math.cos(2 * (lum - lam22)) + 2 * x[0] * math.sin(2 * (lum - lam22)))
    xc22[2] = -5 * A22 * x[2] * (x[2] * x[2] / r / r - 1) * math.cos(lum - lam22)

    return x,xc20,xc30,xc40,xc50,xc60,xc22
#################################################################################

def Atm(x,W):
    Cx = 0.1
    S = 0.944
    a = 1
    A = 0
    hj = 0
    k1 = 0
    k2 = 0
    if (r_r >= 0) and (r_r < 20):
        hj = 0
        A = 1.225e+9
        k1 = 7.825 * 100
        k2 = -263.9e+5

    if (r_r >= 20) and (r_r < 60):
        hj = 20
        A = 0.891e+8
        k1 = 16.37 * 100
        k2 = 44.07e+5

    if (r_r >= 60) and (r_r < 100):
        hj = 60
        A = 2.578e+5
        k1 = 5.905 * 100
        k2 = -256e+5

    if (r_r >= 100) and (r_r < 150):
        hj = 100
        A = 4.061e+2
        k1 = 17.87 * 100
        k2 = 146.9e+5

    if (r_r >= 150) and (r_r < 300):
        hj = 150
        A = 2.13
        k1 = 3.734 * 100
        k2 = 8.004e+5

    if (r_r >= 300) and (r_r < 600):
        hj = 300
        A = 4.764e-2
        k1 = 0.7735 * 100
        k2 = 0.7111e+5

    if (r_r >= 600) and (r_r < 900):
        hj = 600
        A = 8.726e-3
        k1 = 0.928 * 100
        k2 = 0.1831e+5

    if r_r >= 900:
        hj = 900
        A = 6.367 - 4
        k1 = 0.954 * 100
        k2 = 0

    p = A * math.exp(-k1 * (r_r - hj) + k2 * (r_r + hj) * (r_r + hj))
    Cb = Cx * S / 2 / m
    V = math.sqrt(x[3] * x[3] + x[4] * x[4] + x[5] * x[5])
    A = Cb * p * V
    W[0] = -A * (x[3] + W0 * x[1])
    W[1] = -A * (x[4] - W0 * x[0])
    W[2] = -A * x[5]

    return x,W
#################################################################################

def Moon_model(dT, Xmoon):

    T = dT - 2400000.5
    ee = 23.43929111 * math.pi / 180

    L = 218.31617 * math.pi / 180 + 481267.88088 * math.pi / 180 * T - 1.3972 * math.pi / 180 * T
    l = 134.96292 * math.pi / 180 + 477198.867553 * math.pi / 180 * T
    dl = 357.52543 * math.pi / 180 + 35999.04955 * math.pi / 180 * T
    F = 93.27283 * math.pi / 180 + 483202.01873 * math.pi / 180 * T
    D = 297.85027 * math.pi / 180 + 445267.11135 * math.pi / 180 * T

    lam = L + 22640.0 / 3600 * math.pi / 180 * math.sin(l) + 769.0 / 3600 * math.pi / 180 * math.sin(2 * l) - 4586.0 / 3600 * math.pi / 180 * math.sin(l - 2 * D) + 2370.0 / 3600 * math.pi / 180 * math.sin(2 * D) - 668.0 / 3600 * math.pi / 180 * math.sin(dl) - 412.0 / 3600 * math.pi / 180 * math.sin(2 * F) - 212.0 / 3600 * math.pi / 180 * math.sin(2 * l - 2 * D) - 206.0 / 3600 * math.pi / 180 * math.sin(l + dl - 2 * D) + 192.0 / 3600 * math.pi / 180 * math.sin(l + 2 * D) - 165.0 / 3600 * math.pi / 180 * math.sin(dl - 2 * D) + 148. / 3600 * math.pi / 180 * math.sin(dl - l) - 125. / 3600 * math.pi / 180 * math.sin(D) - 110.0 / 3600 * math.pi / 180 * math.sin(l + dl) - 55.0 / 3600 * math.pi / 180 * math.sin(2 * F - 2 * D)
    Bm = 18520. / 3600 * math.pi / 180 * math.sin(F + lam - L + 412.0 / 3600 * math.pi / 180 * math.sin(2 * F) + 541.0 / 3600 * math.pi / 180 * math.sin(dl)) - 526.0 / 3600 * math.pi / 180 * math.sin(F - 2 * D) + 44.0 / 3600 * math.pi / 180 * math.sin(l + F - 2 * D) - 31.0 / 3600 * math.pi / 180 * math.sin(-l + F - 2 * D) - 25.0 / 3600 * math.pi / 180 * math.sin(-2 * l + F) - 23.0 / 3600 * math.pi / 180 * math.sin(dl + F - 2 * D) + 21.0 / 3600 * math.pi / 180 * math.sin(-l + F) + 11.0 / 3600 * math.pi / 180 * math.sin(-dl + F - 2 * D)
    Rm = 385000 - 20905 * math.cos(l) - 3699 * math.cos(2 * D - l) - 2956 * math.cos(2 * D) - 570 * math.cos(2 * l) + 246 * math.cos(2 * l - 2 * D) - 205 * math.cos(dl - 2 * D) - 171 * math.cos(l + 2 * D)
    Xmoon[0] = Rm * math.cos(Bm) * math.cos(lam)
    Xmoon[1] = Rm * (math.cos(Bm) * math.cos(ee) * math.sin(lam) - math.sin(Bm) * math.sin(ee))
    Xmoon[2] = Rm * (math.cos(Bm) * math.sin(ee) * math.sin(lam) - math.sin(Bm) * math.cos(ee))

    return Xmoon
#################################################################################


def Sun_model(Tt,Xsun):

    T = Tt - 2400000.5
    ee = 23.43929111 * math.pi / 180
    L = 218.31617 * math.pi / 180 + 481267.88088 * math.pi / 180 * T - 1.3972 * math.pi / 180 * T
    l = 134.96292 * math.pi / 180 + 477198.867553 * math.pi / 180 * T
    dl = 357.52543 * math.pi / 180 + 35999.04955 * math.pi / 180 * T
    F = 93.27283 * math.pi / 180 + 483202.01873 * math.pi / 180 * T
    D = 297.85027 * math.pi / 180 + 445267.11135 * math.pi / 180 * T
    M = 357.5256 * math.pi / 180 + 35999.049 * math.pi / 180 * T
    lam = 292.9400 * math.pi / 180 + M + 6892.0 / 3600 * math.pi / 180 * math.sin(M) + 72. / 3600 * math.pi / 180 * math.sin(2 * M)
    Rs = (149.619 - 2.499 * math.cos(M) - 0.021 * math.cos(2 * M)) * math.pow(10, 6)
    Xsun[0] = Rs * math.cos(lam)
    Xsun[1] = Rs * math.sin(lam) * math.cos(ee)
    Xsun[2] = Rs * math.sin(lam) * math.sin(ee)

    return Xsun
#################################################################################

def MoonandSunVozm(Xsun, Xmoon,x,G):

    G1 = [0, 0, 0]
    G2 = [0, 0, 0]
    delm = math.sqrt((x[0] - Xmoon[0]) * (x[0] - Xmoon[0]) + (x[1] - Xmoon[1]) * (x[1] - Xmoon[1]) + (x[2] - Xmoon[2]) * (x[2] - Xmoon[2]))
    dels = math.sqrt((x[0] - Xsun[0]) * (x[0] - Xsun[0]) + (x[1] - Xsun[1]) * (x[1] - Xsun[1]) +(x[2] - Xsun[2]) * (x[2] - Xsun[2]))
    Rm = math.sqrt(Xmoon[0] * Xmoon[0] + Xmoon[1] * Xmoon[1] + Xmoon[2] * Xmoon[2])
    Rs = math.sqrt(Xsun[0] * Xsun[0] + Xsun[1] * Xsun[1] + Xsun[2] * Xsun[2])
    G1[0] = MUM * ((Xmoon[0] - x[0]) / math.pow(delm, 3) - Xmoon[0] / math.pow(Rm, 3))
    G1[1] = MUM * ((Xmoon[1] - x[1]) / math.pow(delm, 3) - Xmoon[1] / math.pow(Rm, 3))
    G1[2] = MUM * ((Xmoon[2] - x[2]) / math.pow(delm, 3) - Xmoon[2] / math.pow(Rm, 3))
    G2[0] = MUS * ((Xsun[0] - x[0]) / math.pow(dels, 3) - Xsun[0] / math.pow(Rs, 3))
    G2[1] = MUS * ((Xsun[1] - x[1]) / math.pow(dels, 3) - Xsun[1] / math.pow(Rs, 3))
    G2[2] = MUS * ((Xsun[2] - x[2]) / math.pow(dels, 3) - Xsun[2] / math.pow(Rs, 3))
    G[0] = G1[0] + G2[0]
    G[1] = G1[1] + G2[1]
    G[2] = G1[2] + G2[2]

    return Xsun,Xmoon,x,G
#################################################################################


def Lw(x,Xsun, Gl):

    b = [0,0,0]
    del_ = [0,0,0,]
    xi = [0.0] * 3
    Cp = 0.1
    S = 0.01
    Rz = 6371
    Rs = 696000
    k = 4.56 * math.pow(10, -12)
    A = 1.4959787061 * math.pow(10, 8)
    R = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    Rr = math.sqrt((x[0] - Xsun[0]) * (x[0] - Xsun[0]) + (x[1] - Xsun[1]) * (x[1] - Xsun[1]) + (x[2] - Xsun[2]) * (x[2] - Xsun[2]))
    del_[0] = Xsun[0] - x[0]
    del_[1] = Xsun[1] - x[1]
    del_[2] = Xsun[2] - x[2]
    b[0] = math.acos(-x[0] * del_[0] / R / Rr)
    b[1] = math.acos(-x[1] * del_[1] / R / Rr)
    b[2] = math.acos(-x[2] * del_[2] / R / Rr)
    Ez = math.asin(Rz / R)
    Es = math.asin(Rs / Rr)
    if b[0] >= (Ez + Es):
        xi[0] = 1
    if b[1] >= (Ez + Es):
        xi[1] = 1
    if b[2] >= (Ez + Es):
        xi[2] = 1
    if b[0] <= (Ez - Es):
        xi[0] = 0
    if b[1] <= (Ez - Es):
        xi[1] = 0
    if b[2] <= (Ez - Es):
        xi[2] = 0
    if (b[0] > (Ez - Es)) and (b[0] < (Ez + Es)):
        xi[0] = (b[0] - Ez + Es) / 2 / Es
    if (b[1] > (Ez - Es)) and (b[1] < (Ez + Es)):
        xi[1] = (b[1] - Ez + Es) / 2 / Es
    if (b[2] > (Ez - Es)) and (b[2] < (Ez + Es)):
        xi[2] = (b[2] - Ez + Es) / 2 / Es
        Gl[0] = xi[0] * k * int((A / R) * (A / R)) * Cp * S / m * x[0] / R
    Gl[1] = xi[1] * k * int((A / R) * (A / R)) * Cp * S / m * x[1] / R
    Gl[2] = xi[2] * k * int((A / R) * (A / R)) * Cp * S / m * x[2] / R

    return x,Xsun,Gl
###########################################################

def RP(t,x,dx):

    F = [0,0,0]
    W = [0,0,0]
    Watm = [0,0,0]
    Wms = [0,0,0]
    Wlight = [0,0,0]
    Xmoon = [0,0,0]
    Xsun = [0,0,0]
    Xg = [0,0,0]
    X20 = [0,0,0]
    X30 = [0,0,0]
    X40 = [0,0,0]
    X22 = [0,0,0]
    X50 = [0,0,0]
    X60 = [0,0,0]


    r = math.sqrt(x[0] * x[0] + x[2] * x[2] + x[1] * x[1])
    GESK(x = x, Xg = Xg,lum=lum)
    Moon_model(JD,Xmoon = Xmoon)
    Sun_model(JD, Xsun = Xsun)
    MoonandSunVozm(Xsun = Xsun, Xmoon=Xmoon, x=x,G=Wms)
    Lw(x=x, Xsun=Xsun, Gl=Wlight)
    Atm(x = x, W = Watm)
    C(x=x, Ss= S, xc20=X20, xc30=X30, xc40=X40, xc50=X50, xc60=X60,xc22= X22)
    X22[0] = 0
    X22[1] = 0
    X22[2] = 0

    W[0] = Watm[0] + Wms[0] + Wlight[0] + X20[0] + X30[0] + X40[0] + X22[0] + X50[0] + X60[0]
    W[1] = Watm[1] + Wms[1] + Wlight[1] + X20[1] + X30[1] + X40[1] + X22[1] + X50[1] + X60[1]
    W[2] = Watm[2] + Wms[2] + Wlight[2] + X20[2] + X30[2] + X40[2] + X22[2] + X50[2] + X60[2]
    dx[0] = x[3]
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = -MU * x[0] / math.pow(r, 3) + W[0]
    dx[4] = -MU * x[1] / math.pow(r, 3) + W[1]
    dx[5] = -MU * x[2] / math.pow(r, 3) + W[2]

    return x, dx
#################################################################################

def Runge(t, x, xk, h):

    k1 = [0] * 6
    k1t= [0] * 6
    k2= [0] * 6
    k2t= [0] * 6
    k3= [0] * 6
    k3t= [0] * 6
    k4= [0] * 6
    RP(t = t, x=x, dx=k1)
    for i in range(0,6,1) :
        k1t[i] = x[i] + (h / 2) * k1[i]
    RP(t + h / 2, x=k1t, dx=k2)
    for i in range(0,6,1) :
        k2t[i] = x[i] + (h / 2) * k2[i]
    RP(t + h / 2, x = k2t, dx = k3)
    for i in range(0,6,1) :
        k3t[i] = x[i] + h * k3[i]
    RP(t + h, x = k3t, dx = k4)
    for i in range(0,6,1) :
        xk[i] = x[i] + (h / 6) * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i])

    return x,xk
#################################################################################

def el(x, x1):

    C = [0,0,0]
    Ff = [0,0,0]
    C[0] = x[1] * x[5] - x[2] * x[4]
    C[1] = x[2] * x[3] - x[0] * x[5]
    C[2] = x[0] * x[4] - x[1] * x[3]

    RRR = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    Ff[0] = C[2] * x[4] - C[1] * x[5] - MU * x[0] / RRR
    Ff[1] = C[0] * x[5] - C[2] * x[3] - MU * x[1] / RRR
    Ff[2] = C[1] * x[3] - C[0] * x[4] - MU * x[2] / RRR
    CC = math.sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2])
    FF = math.sqrt(Ff[0] * Ff[0] + Ff[1] * Ff[1] + Ff[2] * Ff[2])
    p = CC * CC / MU
    x1[1] = FF / MU
    x1[0] = p / (1 + x1[1])

    x1[2] = math.acos(C[2] / CC)

    if x1[2] != 0:
        Om_ = math.acos(-C[1] / (CC * math.sin(x1[2])))
        if (C[0] * math.sin(x1[2]) / CC) < 0:
            x1[5] = 2 * math.pi - Om_
        else :
            x1[5] = Om_
        u_ = math.acos((x[0] * math.cos(x1[5]) + x[1] * math.sin(x1[5])) / RRR)
        if x[2] < 0:
            x1[3] = 2 * math.pi - u_
        else :
            x1[3] = u_

    if x1[2] == 0:
        u_ = math.acos(x[0] / RRR)
        if x[1] >= 0:
            x1[3] = u_
        else:
            x1[3] = 2 * math.pi - u_
        x1[5] = x1[3]

    if (x1[1] > 0) and (x1[1] < 1):
        if math.cos(x1[2]) == 0:
            Sw = Ff[2] / FF
            Cw = Ff[1] / FF
        else :
            Sw = (-Ff[0] / FF * math.sin(x1[5]) + Ff[1] / FF * math.cos(x1[5])) / math.cos(x1[2])
            Cw = Ff[0] / FF * math.cos(x1[5]) + Ff[1] / FF * math.sin(x1[5])

        if Sw >= 0:
            x1[4] = math.acos(Cw)
        else :
            x1[4] = 2 * math.pi - math.acos(Cw)

    if x1[1] == 0:
        x1[4] = x1[3]

    return x,x1
#################################################################################

def Startdata(x, Rp, e, i, u, w, U):

    x1[0] = Rp
    print(x[0])
    x1[1] = e
    x1[2] = i * math.pi / 180
    x1[3] = u * math.pi / 180
    x1[4] = w * math.pi / 180
    x1[5] = U * math.pi / 180

    return x
#################################################################################

def GetXYZ(x,x1):

    p = x1[0] * (1 + x1[1])
    v = x1[3] - x1[4]
    r = p / (1 + x1[1] * math.cos(v))
    dr = math.sqrt(MU / p) * x1[1] * math.sin(v)
    dv = math.sqrt(MU * p) / r / r

    x[0] = r * (math.cos(x1[3]) * math.cos(x1[5]) - math.sin(x1[3]) * math.sin(x1[5]) * math.cos(x1[2]))
    x[1] = r * (math.cos(x1[3]) * math.sin(x1[5]) + math.sin(x1[3]) * math.cos(x1[5]) * math.cos(x1[2]))
    x[2] = r * math.sin(x1[3]) * math.sin(x1[2])
    x[3] = dr * x[0] / r - r * (math.sin(x1[3]) * math.cos(x1[5]) + math.cos(x1[3]) * math.sin(x1[5]) * math.cos(x1[2])) * dv
    x[4] = dr * x[1] / r - r * (math.sin(x1[3]) * math.sin(x1[5]) - math.cos(x1[3]) * math.cos(x1[5]) * math.cos(x1[2])) * dv
    x[5] = dr * x[2] / r + r * dv * math.cos(x1[3]) * math.sin(x1[2])

    return x,x1
#################################################################################


def main():

    xk = [0,0,0,0,0,0,0]
    x= [0,0,0,0,0,0,0]
    t = 0
    rp = [None]
    global x1
    Startdata(x = x1, Rp=7852.33737599277, e = 0.0022291, i = 82.5088, u = 0.107127599236494,w =  141.9354, U= 150.9592)
    global m
    m = 280
    tk = 2592000
    h = 60
    p = x1[0] * (1 + x1[1])
    v = x1[3] - x1[4]
    r = p / (1 + x1[1] * math.cos(v))
    global r_r
    r_r = r

    GetXYZ(x=x, x1 = x1)

    print("t Rp e i u w U\n")
    sh = 0
    rp[0] = x1[0]
    times = []
    Rpgraph = []
    egraph = []
    igraph = []
    ugraph = []
    wgraph = []
    Ugraph = []
    while t <= tk:
        #Вывод параметров в консоль
        print( t, x1[0], x1[1], x1[2] / math.pi * 180, x1[3] / math.pi * 180,x1[4] / math.pi * 180, x1[5] / math.pi * 180)
        if t % 60000==0:
            times.append(t)
            Rpgraph.append(x1[0])
            egraph.append(x1[1])
            igraph.append(x1[2] / math.pi * 180)
            ugraph.append(x1[3] / math.pi * 180)
            wgraph.append(x1[4] / math.pi * 180)
            Ugraph.append(x1[5] / math.pi * 180)

        JulyData(t = t, JD=JD)
        TimeStar(t, JD)
        Runge(t = t, x=x, xk=xk, h=h)
        for i in range(0,6,1) :
            x[i] = xk[i]
        el(x=x, x1=x1)
        sh += 1
        t += h
    #Вывод графиков
    plt.xlabel("t")  # ось абсцисс
    plt.ylabel("Rp")  # ось ординат
    plt.grid()  # включение отображение сетки
    plt.plot(times, Rpgraph)  # построение графика
    plt.show()
    plt.xlabel("t")  # ось абсцисс
    plt.ylabel("e")  # ось ординат
    plt.grid()  # включение отображение сетки
    plt.plot(times, egraph)  # построение графика
    plt.show()
    plt.xlabel("t")  # ось абсцисс
    plt.ylabel("i")  # ось ординат
    plt.grid()  # включение отображение сетки
    plt.plot(times, igraph)  # построение графика
    plt.show()
    plt.xlabel("t")  # ось абсцисс
    plt.ylabel("u")  # ось ординат
    plt.grid()  # включение отображение сетки
    plt.plot(times, ugraph)  # построение графика
    plt.show()
    plt.xlabel("t")  # ось абсцисс
    plt.ylabel("w")  # ось ординат
    plt.grid()  # включение отображение сетки
    plt.plot(times, wgraph)  # построение графика
    plt.show()
    plt.xlabel("t")  # ось абсцисс
    plt.ylabel("U")  # ось ординат
    plt.grid()  # включение отображение сетки
    plt.plot(times, Ugraph)  # построение графика
    plt.show()
    input()



main()



