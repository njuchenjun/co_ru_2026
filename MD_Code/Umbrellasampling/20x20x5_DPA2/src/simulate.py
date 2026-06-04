import src.state as state
import src.md as md

def run_simulation():
    state.uout.write('Umbrella Sampling starts...')
    state.uout.write(f"Sampling Window : {state.us_ius:>6}\n")
    state.uout.write(f"Equilibrium position : {state.us_z0:>12.5f}, Force constant : {state.us_fc:>12.5f}, eV/A/A\n")
    state.uout.write(f"{'Step:':>12} {'Temperature (K)':>20} {'Potential (eV)':>20} {'Biased V (eV)':>20} {'Total energy (eV)':>20}\n")
    state.uout.write("-" * 92 + "\n")

    avg_T = 0.0
    var_T = 0.0
    
    # 生成文件名：US-aXXX 和 tra-aXXX，其中 XXX 为窗口编号（右对齐，前补零）
    winname = f"US-a{state.us_ius:03d}"
    # tracename = f"tra-a{state.us_ius:03d}"

    # 打开文件
    state.uwin = open(winname, 'w',buffering=1)
    # stuwc = open(tracename, 'w')

    for i_step in range(1, state.nstep + 1):
        md.MD_VV(i_step,state.zmc,state.us_z0)

        avg_T += state.tempnow
        var_T += (state.tempnow - state.Temperature) ** 2

        if i_step % state.outstep == 0:
            pot_eV = state.Ev / 100.0 / 96.4853
            tot_eV = (state.Ev + state.Ek) / 100.0 / 96.4853
            state.uout.write(f"{i_step:>12} {state.tempnow:>20.6f} {pot_eV:>20.8e} {state.Ebias:>20.8e} {tot_eV:>20.8e}\n")
            state.uwin.write(f"{i_step-1:>15} {state.zmc:>16.8f}\n")

        if i_step % (state.outstep) == 0:
            #state.utrjtxt.write(f"{state.natom-(state.supercell[0]*state.supercell[1]*state.fixed)}\n")
            state.utrjtxt.write(f"{2}\n")
            state.utrjtxt.write(f"{i_step:>2} {pot_eV:>18.8e}\n")
            for i in range((state.supercell[0]*state.supercell[1]*5),state.supercell[0]*state.supercell[1]*5+2):
                state.utrjtxt.write(f"{state.element[i]:>2}  {state.Cmat[0,i]:>16.6e} {state.Cmat[1,i]:>16.6e} {state.Cmat[2,i]:>16.6e}\n")

    state.uout.write("-" * 92 + "\n")
    state.uout.write("End of simulation.\n\n")

    avg_T /= state.nstep
    var_T /= state.nstep
    expt_var_T = 2.0 * state.Temperature ** 2 / state.nfr

    state.uout.write(f"Expected temperature: {state.Temperature:>14.3f}\n")
    state.uout.write(f"Averaged temperature: {avg_T:>14.3f}\n\n")
    state.uout.write(f"Expected variance of temperature: {expt_var_T:>14.3f}\n")
    state.uout.write(f"Variance of temperature: {var_T:>14.3f}\n\n")
    
    state.uwin.close()
