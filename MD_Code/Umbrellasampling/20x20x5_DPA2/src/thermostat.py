# src/thermostat.py
import math
import src.state as state
import src.const as const
from src.utils import gaussian_knuth

NNOS = 6
NYOSH = 5
NRESN = 2

# 恒温器数组
WDTI2 = [0.0] * NYOSH
WDTI4 = [0.0] * NYOSH
WDTI8 = [0.0] * NYOSH
Qmass = [0.0] * NNOS
Glogs = [0.0] * NNOS
Vlogs = [0.0] * NNOS
Xlogs = [0.0] * NNOS
GNKT = 0.0
GKT = 0.0

def iniNHC():
    global GNKT, GKT, Qmass, Vlogs, Xlogs, Glogs, WDTI2, WDTI4, WDTI8
    GNKT = float(state.nfr) * state.Temperature * const.ektoint
    GKT = state.Temperature * const.ektoint
    Nfreq = const.fcmtoint * state.fcpl * 2.0 * const.pi

    Qmass[0] = GNKT / (Nfreq * Nfreq)
    for i in range(1, NNOS):
        Qmass[i] = GKT / (Nfreq * Nfreq)

    # Yoshida 系数
    tmp = [0.0] * NYOSH
    for i in range(NYOSH):
        tmp[i] = 1.0 / (4.0 - 1.5874010519682)
    tmp[2] = 1.0 - 4.0 * tmp[0]

    for i in range(NYOSH):
        WDTI2[i] = state.timestep * tmp[i] / (2.0 * NRESN)
        WDTI4[i] = state.timestep * tmp[i] / (4.0 * NRESN)
        WDTI8[i] = state.timestep * tmp[i] / (8.0 * NRESN)

    if state.debug == 1:
        state.udeb.write("Initialize NHC thermostat\n")
        state.udeb.write(f"GNKT={GNKT}\n")
        state.udeb.write(f"GKT={GKT}\n")
        state.udeb.write("WDTI2\n")
        state.udeb.write(" ".join(str(WDTI2[i]) for i in range(NYOSH)) + "\n")
        state.udeb.write("WDTI4\n")
        state.udeb.write(" ".join(str(WDTI4[i]) for i in range(NYOSH)) + "\n")
        state.udeb.write("WDTI8\n")
        state.udeb.write(" ".join(str(WDTI8[i]) for i in range(NYOSH)) + "\n")

    Glogs[0] = 0.0
    Xlogs[0] = 1.0
    for i in range(1, NNOS):
        Glogs[i] = 0.0
        Xlogs[i] = 1.0

    Vlogs[0] = math.sqrt(2.0 * GNKT / Qmass[0])
    for i in range(1, NNOS):
        Vlogs[i] = gaussian_knuth() / math.sqrt(Qmass[i])

    ttt = sum(Vlogs[i]**2 * Qmass[i] for i in range(1, NNOS))
    if ttt > 1e-9:
        factor = math.sqrt(GKT * (NNOS - 1) / ttt)
        for i in range(1, NNOS):
            Vlogs[i] *= factor

def NHCINT():
    """积分 Nose–Hoover 链从 t=0 到 t=timestep/2"""
    global GNKT, GKT, Qmass, Glogs, Vlogs, Xlogs, WDTI2, WDTI4, WDTI8
    # 延迟导入以避免循环依赖
    from src.md import Calctemp, RescaleV

    # 计算当前动能
    Calctemp()
    AKIN = 2.0 * state.Ek

    # 调试输出
    if state.debug == 1:
        state.udeb.write("NHCINT\n")
        state.udeb.write("GLOGS\n")
        state.udeb.write(" ".join(str(Glogs[j]) for j in range(NNOS)) + "\n")
        state.udeb.write("QMASS\n")
        state.udeb.write(" ".join(str(Qmass[j]) for j in range(NNOS)) + "\n")
        state.udeb.write("XLOGS\n")
        state.udeb.write(" ".join(str(Xlogs[j]) for j in range(NNOS)) + "\n")
        state.udeb.write("VLOGS\n")
        state.udeb.write(" ".join(str(Vlogs[j]) for j in range(NNOS)) + "\n")

    cscale = 1.0

    for IRESN in range(NRESN):
        for IYOSH in range(NYOSH):
            # ---------- 更新恒温器力 ----------
            # 最后一个链节
            Glogs[NNOS-1] = (Qmass[NNOS-2] * Vlogs[NNOS-2]**2 - GKT) / Qmass[NNOS-1]

            # 中间链节（从 NNOS-2 到 1）
            for INOS in range(NNOS-2, 0, -1):
                i = INOS - 1   # 0-based 索引
                Glogs[i] = (Qmass[i-1] * Vlogs[i-1]**2 - GKT) / Qmass[i]

            # 第一个链节
            Glogs[0] = (AKIN - GNKT) / Qmass[0]

            # ---------- 更新恒温器速度（第一部分） ----------
            # 最后一个链节
            Vlogs[NNOS-1] += Glogs[NNOS-1] * WDTI4[IYOSH]

            # 其余链节（从 NNOS-1 向下到 1）
            for INOS in range(1, NNOS):
                idx = NNOS - 1 - INOS   # 对应 Fortran 的 NNOS-INOS
                AA = math.exp(-WDTI8[IYOSH] * Vlogs[idx+1])
                Vlogs[idx] = Vlogs[idx] * AA * AA + WDTI4[IYOSH] * Glogs[idx] * AA

            # ---------- 更新粒子速度和动能 ----------
            AA = math.exp(-WDTI2[IYOSH] * Vlogs[0])
            cscale *= AA
            AKIN *= AA * AA

            if state.debug == 1:
                state.udeb.write(f"AA={AA}\n")

            # ---------- 更新恒温器位置 ----------
            for INOS in range(NNOS):
                Xlogs[INOS] += Vlogs[INOS] * WDTI2[IYOSH]

            # ---------- 再次更新恒温器力 ----------
            Glogs[0] = (AKIN - GNKT) / Qmass[0]

            # ---------- 更新恒温器速度和力（第二部分） ----------
            for INOS in range(NNOS - 1):
                AA = math.exp(-WDTI8[IYOSH] * Vlogs[INOS+1])
                Vlogs[INOS] = Vlogs[INOS] * AA * AA + WDTI4[IYOSH] * Glogs[INOS] * AA
                Glogs[INOS+1] = (Qmass[INOS] * Vlogs[INOS]**2 - GKT) / Qmass[INOS+1]

            # 最后一个链节速度更新（第二次）
            Vlogs[NNOS-1] += Glogs[NNOS-1] * WDTI4[IYOSH]

            if state.debug == 1:
                state.udeb.write(f"Scale={cscale}\n")

    if state.debug == 1:
        state.udeb.write(f"NHCINT Rescale factor {cscale}\n")

    # 应用速度缩放因子
    RescaleV(cscale)


def NHC_VV(job):
    from src.md import MD_VV_half
    if job == 1:
        NHCINT()
        MD_VV_half()
    elif job == 2:
        MD_VV_half()
        NHCINT()