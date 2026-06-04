#!/usr/bin/env python3
# pes.py - 势能面接口，使用 DeePMD-kit 模型

import numpy as np
from deepmd.infer import DeepPot
import src.state as state

# 全局模型对象（单例模式，只加载一次）
_deepmd_model = None
_MODEL_PATH = "models/dpa2_m_OC20M247f.pth"  # 根据你的模型路径修改

def _get_model():
    """懒加载 DeePMD 模型"""
    global _deepmd_model
    if _deepmd_model is None:
        print(f"Loading DeePMD model from {_MODEL_PATH}...")
        _deepmd_model = DeepPot(_MODEL_PATH)
        print("Model loaded successfully.")
    return _deepmd_model

def pes_nCO_RuHD(natom, Cmat, egy, dedx):
    """
    使用 DeePMD-kit 模型计算势能和力

    输入:
        natom: int, 原子数
        Cmat: (3, natom) numpy array, 原子坐标（单位 Å）
        egy: 长度为1的数组，用于存储势能（eV）
        dedx: (3, natom) numpy array，用于存储力（负梯度，单位 eV/Å）
    """
    # 获取模型
    dp = _get_model()

    # 1. 构建坐标数组：shape (1, natom, 3)
    #    Cmat 当前是 (3, natom)，需要转置为 (natom, 3)，然后增加批次维度
    coord = Cmat.T.reshape(1, natom, 3)

    # 2. 构建晶胞（盒子）数组：shape (1, 3:, 3)
    #    使用参考中的构造方式，从模拟参数获取晶格常数和超胞尺寸
    #    假设 x 和 y 方向的超胞尺寸分别为 supercell[0] 和 supercell[1]
    nx = state.supercell[0]
    ny = state.supercell[1]
    a = state.lconst          # 晶格常数（单位 Å）
    
    # 二维六角格子基矢（未乘超胞）
    # 参考： lattice_vector(1,1) = sqrt(3)/2, lattice_vector(2,1) = -0.5
    #        lattice_vector(1,2) = 0.0,   lattice_vector(2,2) = 1.0
    # 乘以超胞尺寸后，盒子的 x 方向基矢为 (sqrt(3)/2 * a * nx, -0.5 * a * nx, 0)
    # y 方向基矢为 (0, a * ny, 0)
    # z 方向使用模型训练时的固定值（根据参考，应为 23.5301303864）
    # 如果你有 z 方向长度，也可以从 state.zbox 读取
    zbox = 23.5301303864   # 可根据实际模型修改
    cell = np.array([
        [np.sqrt(3.0)/2.0 * a * nx, -0.5 * a * nx, 0.0],
        [0.0, a * ny, 0.0],
        [0.0, 0.0, zbox]
    ]).reshape(1, 3, 3)

    # 3. 构建原子类型列表 atype
    #    假设 DeePMD 模型训练时的 type_map 为 ["Ru", "C", "O"]
    #    即 Ru = 0, C = 1, O = 2
    #    前 48 个原子是表面 Ru，后面的原子为 CO 分子对（C 和 O 交替）
    #    根据元素符号映射
    type_map = {"Ru": 0, "C": 1, "O": 2}
    atype = [type_map[elem] for elem in state.element]  # 从全局 state.element 获取

    # 4. 调用 DeePMD 模型
    #    返回：e (energy), f (force), v (virial)
    e, f, v = dp.eval(coord, cell, atype)

    # 5. 将结果写入输出数组
    #    能量：e 的形状可能是 (1,)，直接取标量
    egy[0] = e[0] if isinstance(e, (list, np.ndarray)) and len(e) == 1 else e

    #    力：f 的形状是 (1, natom, 3)，需要转置为 (3, natom)
    force = f[0].T   # shape (3, natom)
    dedx[:, :] = force[:, :]

    for i in range(0, (nx * ny * state.fixed)):
        dedx[0,i] = 0.0
        dedx[1,i] = 0.0
        dedx[2,i] = 0.0

    # 可选：如果你需要记录 virial，可以存储在全局状态中
    # state.virial = v[0]  # shape (9,)
