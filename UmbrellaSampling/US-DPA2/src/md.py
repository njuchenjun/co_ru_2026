import math
import numpy as np
import src.state as state
import src.const as const
import src.thermostat as thermostat
import src.pes as pes

def MD_VV_half():
    hstep = 0.5 * state.timestep
    for i in range(state.supercell[0]*state.supercell[1]*state.fixed,state.natom):
        for j in range(3):
            state.Vmat[j, i] += state.Fmat[j, i] / state.Mass[i] * hstep

def Calctemp():
    state.Ek = 0.0
    for i in range(state.supercell[0]*state.supercell[1]*state.fixed,state.natom):
        for j in range(3):
            state.Ek += state.Mass[i] * state.Vmat[j, i] ** 2
    state.tempnow = state.Ek / const.ektoint / float(state.nfr)
    state.Ek *= 0.5

def RescaleV(cscale):
    for i in range(state.supercell[0]*state.supercell[1]*state.fixed,state.natom):
        for j in range(3):
            state.Vmat[j, i] *= cscale

def MD_VV(i_step,zmc,zeq):
    # 第一半速度更新
    if state.entype == 1:
        thermostat.NHC_VV(1)
    elif state.entype == 0:
        MD_VV_half()

    # 坐标更新
    for i in range(state.supercell[0]*state.supercell[1]*state.fixed,state.natom):
        for j in range(3):
            state.Cmat[j, i] += state.Vmat[j, i] * state.timestep

    # 势能面调用
    egy = np.array([0.0])
    dedx = np.zeros((3, state.natom))
    pes.pes_nCO_RuHD(state.natom, state.Cmat, egy, dedx)
    
    # 计算 CO 质心坐标（集体变量）
    # 使用质量比或直接计算
    state.zmc = (state.Cmat[2, state.us_iC] * state.Mass[state.us_iC] +
           state.Cmat[2, state.us_iO] * state.Mass[state.us_iO]) / (state.Mass[state.us_iC] + state.Mass[state.us_iO])
    # 或者使用预计算的质量比：
    # zmc = state.Cmat[2, state.us_iC] * state.us_rmC + state.Cmat[2, state.us_iO] * state.us_rmO

    # 计算偏置势能和力
    dz = state.zmc - state.us_z0
    state.Ebias = 0.5 * state.us_fc * dz * dz
    egy[0] += state.Ebias

    # 添加偏置力到 dedx（负梯度方向）
    dedx[2, state.us_iC] -= state.us_fc * dz * state.us_rmC
    dedx[2, state.us_iO] -= state.us_fc * dz * state.us_rmO

    state.Ev = egy[0] * const.eeVtoint
    for i in range(state.supercell[0]*state.supercell[1]*state.fixed,state.natom):
        for j in range(3):
            state.Fmat[j, i] = dedx[j, i] * const.eeVtoint
   # print(f"Step {i_step}: Fmat max = {np.max(np.abs(state.Fmat)):.6e}")

    # 第二半速度更新
    if state.entype == 1:
        thermostat.NHC_VV(2)
    elif state.entype == 0:
        MD_VV_half()

    Calctemp()

    # CO 反弹处理（示例）
 #   num_CO = (state.natom - 48) // 2
 #   for i in range(num_CO):
 #       idx_C = 48 + i * 2
 #       idx_O = 48 + i * 2 + 1
 #       if state.Cmat[2, idx_C] > state.zupper and state.Vmat[2, idx_C] > 0.0:
 #           state.Vmat[2, idx_C] = -state.Vmat[2, idx_C]
 #           state.Vmat[2, idx_O] = -state.Vmat[2, idx_O]
 #           state.uout.write(f"Desorption (ps) / i-th CO {i_step * state.timestep} {i+1}\n")
