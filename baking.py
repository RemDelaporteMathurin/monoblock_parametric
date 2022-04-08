from sim_monoblock import simulate_monoblock_baking
import numpy as np

if __name__ == "__main__":
    T_baking = np.linspace(340, 620, num=6, endpoint=True)
    for T in T_baking:
        simulate_monoblock_baking(
            T_baking=T,
            folder="solution_baking/noninstant_recomb_on_coolant/T={:.2e}".format(T))
