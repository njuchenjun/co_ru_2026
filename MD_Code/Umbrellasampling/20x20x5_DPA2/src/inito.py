import math
import random
import numpy as np
import src.state as state
import src.const as const
import src.utils as utils
import src.thermostat as thermostat
import src.md as md

def read_input():
    with open('INPUT', 'r') as f:
        # 读取所有非空行，但保留行顺序
        lines = [line.strip() for line in f if line.strip() != '']

    # 使用索引遍历
    idx = 0
    # 1. Ensemble
    if lines[idx] == 'Ensemble':
        idx += 1
        state.ensemble = lines[idx]
    else:
        raise ValueError("Expected 'Ensemble' line")
    idx += 1

    # 2. Number of steps
    if lines[idx] == 'Number of steps':
        idx += 1
        state.nstep = int(lines[idx])
    else:
        raise ValueError("Expected 'Number of steps' line")
    idx += 1

    # 3. Output steps
    if lines[idx] == 'Output steps':
        idx += 1
        state.outstep = int(lines[idx])
    else:
        raise ValueError("Expected 'Output steps' line")
    idx += 1

    # 4. Time step(ps)
    if lines[idx] == 'Time step(ps)':
        idx += 1
        state.timestep = float(lines[idx])
    else:
        raise ValueError("Expected 'Time step(ps)' line")
    idx += 1

    # 5. Temperature(K)
    if lines[idx] == 'Temperature(K)':
        idx += 1
        state.Temperature = float(lines[idx])
    else:
        raise ValueError("Expected 'Temperature(K)' line")
    idx += 1

    # 6. Couple frequency(cm-1)
    if lines[idx] == 'Couple frequency(cm-1)':
        idx += 1
        state.fcpl = float(lines[idx])
    else:
        raise ValueError("Expected 'Couple frequency(cm-1)' line")
    idx += 1

    # 7. Lattice constant
    if lines[idx] == 'Lattice constant':
        idx += 1
        state.lconst = float(lines[idx])
    else:
        raise ValueError("Expected 'Lattice constant' line")
    idx += 1

    # 8. Super Cell
    if lines[idx] == 'Super Cell':
        idx += 1
        # 一行包含两个整数
        parts = lines[idx].split()
        if len(parts) != 2:
            raise ValueError("Super Cell line must contain two integers")
        state.supercell[0] = int(parts[0])
        state.supercell[1] = int(parts[1])
    else:
        raise ValueError("Expected 'Super Cell' line")
    idx += 1

    # 9. Number of atoms (adsorption)
    if lines[idx] == 'Number of atoms (adsorption)':
        idx += 1
        state.natom = int(lines[idx])
    else:
        raise ValueError("Expected 'Number of atoms (adsorption)' line")
    idx += 1

    # 分配数组
    state.element = [''] * state.natom
    state.Mass = [0.0] * state.natom
    state.Cmat = np.zeros((3, state.natom), dtype=float)
    state.Cmatsingle = np.zeros((3, state.natom), dtype=float)
    state.Vmat = np.zeros((3, state.natom), dtype=float)
    state.Fmat = np.zeros((3, state.natom), dtype=float)

    # 10. 原子数据标题行
    # 跳过标题行 "Elem   Mass       X         Y        Z"
    if lines[idx] == 'Elem   Mass       X         Y        Z':
        idx += 1
    else:
        raise ValueError("Expected atomic data header line")

    # 读取原子数据
    for i in range(state.natom):
        parts = lines[idx].split()
        if len(parts) < 5:
            raise ValueError(f"Expected 5 fields for atom {i+1}")
        state.element[i] = parts[0]
        state.Mass[i] = float(parts[1])
        state.Cmat[0, i] = float(parts[2])
        state.Cmat[1, i] = float(parts[3])
        state.Cmat[2, i] = float(parts[4])
        idx += 1

    # 11. Upper bound of z-direction
    if lines[idx] == 'Upper bound of z-direction':
        idx += 1
        state.zupper = float(lines[idx])
    else:
        raise ValueError("Expected 'Upper bound of z-direction' line")
    idx += 1

    # 12. Lower bound of potential energy(eV)
    if lines[idx] == 'Lower bound of potential energy(eV)':   # 注意原文件可能是拼写错误 "Louer" 还是 "Lower"？
        # 这里使用用户提供的 "Louer" 保持一致，但也可以兼容两种
        idx += 1
        state.Vlower = float(lines[idx])
    else:
        raise ValueError("Expected 'Lower bound of potential energy(eV)' line")
    idx += 1
    
    if idx < len(lines) and lines[idx] == 'Umbrella sampling':
        idx += 1

    # 读取 Window index
    if idx < len(lines) and lines[idx] == 'Window index':
        idx += 1
        state.us_ius = int(lines[idx])
    else:
        raise ValueError("Expected 'Window index' line")
    idx += 1

    # 读取 Equilibrium z0 value (A)
    if idx < len(lines) and lines[idx] == 'Equilibrium z0 value (A)':
        idx += 1
        state.us_z0 = float(lines[idx])
    else:
        raise ValueError("Expected 'Equilibrium z0 value (A)' line")
    idx += 1

    # 读取 Force constant (单位可能是 kcal/mol/Å²，需要转换为 eV/Å²)
    # 标题行可能包含单位说明，我们匹配其开头 "Force constant"
    if idx < len(lines) and lines[idx].startswith('Force constant'):
        idx += 1
        fc_kcal = float(lines[idx])
        # 转换：1 kcal/mol = 0.043364 eV
        state.us_fc = fc_kcal / const.eeVtokcal
    else:
        raise ValueError("Expected 'Force constant' line")
    idx += 1
        
def setup_output():
    state.uout = open('traj.log', 'w',buffering=1)
    state.utrjtxt = open('traj.txt', 'w',buffering=1)
    if state.debug == 1:
        state.udeb = open('DEBUG', 'w',buffering=1)

def write_header():
    state.uout.write("  MD simulation for CO  Umbrella Sampling\n\n")
    state.uout.write(f"  Ensemble: {state.ensemble}\n\n")
    state.uout.write(f"  Number of steps: {state.nstep:14d}\n\n")
    state.uout.write(f"  Output interval: {state.outstep:14d}\n\n")
    state.uout.write(f"  Time step: {state.timestep:14.6e}\n\n")
    state.uout.write(f"  Temperature: {state.Temperature:12.3f}\n\n")
    state.uout.write(f"  Lattice constant: {state.lconst:14.6e}\n\n")
    state.uout.write(f"  Supercell: {state.supercell[0]:8d} {state.supercell[1]:8d}\n\n")
    state.uout.write(f"  Number of atoms {state.natom:14d}\n\n")
    state.uout.write("  Coordinates\n")
    for i in range(state.natom):
        state.uout.write(f"  {state.element[i]:2s}  {state.Cmat[0,i]:14.6e} {state.Cmat[1,i]:14.6e} {state.Cmat[2,i]:14.6e}\n")
    state.uout.write("\n")
    state.uout.write(f"  Upper bound of z direction: {state.zupper:16.8f}\n\n")
    state.uout.write(f"  Lower bound of potential energy: {state.Vlower:16.8f} eV\n")
    state.Vlower = state.Vlower * const.eeVtoint
    state.uout.write(f"  Lower bound of potential energy: {state.Vlower:18.8e}\n\n")

    state.nfr = (state.natom -(state.supercell[0]*state.supercell[1]*state.fixed))* 3
    state.uout.write(f"  Number of freedom: {state.nfr:8d}\n\n")
    if state.debug == 0:
        state.uout.write("  Debug:          disable\n\n")
    else:
        state.uout.write("  Debug:           enable\n\n")

def setup_lattice():
    state.lattice_vector[0, 0] = math.sqrt(3.0) / 2.0
    state.lattice_vector[1, 0] = -0.5
    state.lattice_vector[0, 1] = 0.0
    state.lattice_vector[1, 1] = 1.0
    state.lsc[0] = float(state.supercell[0]) * state.lconst
    state.lsc[1] = float(state.supercell[1]) * state.lconst

def init_velocities():
    for i in range((state.supercell[0]*state.supercell[1]*state.fixed),state.natom):
        for j in range(3):
            state.Vmat[j, i] = utils.fast_gaussian() / math.sqrt(state.Mass[i])
    md.Calctemp()
    cscale = math.sqrt(state.Temperature / state.tempnow)
    md.RescaleV(cscale)

def initialize():
    read_input()
    
    state.us_iC = 2000
    state.us_iO = 2001
    mass_C = state.Mass[state.us_iC]
    mass_O = state.Mass[state.us_iO]
    total_mass = mass_C + mass_O
    state.us_rmC = mass_C / total_mass
    state.us_rmO = mass_O / total_mass
    print(total_mass)
    
    setup_output()
    if state.ensemble == 'NVE':
        state.entype = 0
    elif state.ensemble == 'NVT':
        state.entype = 1
    write_header()
    setup_lattice()
    init_velocities()
    thermostat.iniNHC()
    

    

