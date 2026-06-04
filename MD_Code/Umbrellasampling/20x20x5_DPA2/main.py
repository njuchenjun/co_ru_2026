#!/usr/bin/env python3
from src.init import initialize
from src.simulate import run_simulation
import src.state as state

if __name__ == "__main__":
    # 可选：设置随机种子
    # import random; random.seed(12345)
    initialize()
    run_simulation()
    # 关闭文件
    state.uout.close()
    state.utrjtxt.close()
    if state.debug == 1 and state.udeb is not None:
        state.udeb.close()