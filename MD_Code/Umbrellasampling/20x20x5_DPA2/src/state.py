import numpy as np

# 模拟参数
ensemble = None          # 'NVE' 或 'NVT'
nstep = None
outstep = None
timestep = None
Temperature = None
fcpl = None
lconst = None
supercell = [0, 0]
natom = None
zupper = None
Vlower = 0.0
debug = 0

# 原子数据
element = []             # 元素符号列表
Mass = []                # 原子质量列表
Cmat = None              # 坐标 (3, natom)
Cmatsingle = None        # 临时坐标副本
Vmat = None              # 速度 (3, natom)
Fmat = None              # 力 (3, natom)
fixed = 3 

# 能量与温度
Ev = 0.0
Ek = 0.0
Ebias = 0.0
tempnow = 0.0
nfr = 0

# 晶格
lattice_vector = np.zeros((2, 2), dtype=float)
lsc = [0.0, 0.0]

# 系综类型
entype = 0               # 0=NVE, 1=NVT

# 文件对象
uinp = None
uout = None
udeb = None
utrjtxt = None
uwin = None

# ==================== 伞形采样相关变量 ====================
us_ius = 0          # 窗口编号（0表示不进行伞形采样）
us_z0 = 0.0         # 平衡位置 (Å)
us_fc = 0.0         # 力常数 (eV/Å²)
us_atom_index = -1  # 集体变量对应的原子索引（0-based，-1表示未指定）

us_iC = -1          # CO 分子中 C 原子的索引
us_iO = -1          # CO 分子中 O 原子的索引
us_rmC = 0.0        # C 原子质量比
us_rmO = 0.0        # O 原子质量比

zmc = 0.0
