from math import cos, sin

DEG_RAD = 0.017453293

def dircos(incl, decl, azim):
    incl_x = incl * DEG_RAD
    decl_x = decl * DEG_RAD
    azim_x = azim * DEG_RAD

    a = cos(incl_x) * cos(decl_x - azim_x)
    b = cos(incl_x) * sin(decl_x - azim_x)
    c = sin(incl_x)

    return a, b, c