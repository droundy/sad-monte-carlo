#!/usr/bin/python3

import numpy as np
from subprocess import run
import system

run(['cargo', 'build', '--release', '--bin',
     'replicas', '--bin', 'histogram', '--bin', 'tempering'], check=True)

max_iter_default = 1e12


def rq(name, cmd, cpus):
    run(f'rq run -c {cpus} --max-output=30 -R -J'.split() +
        [name, '--']+cmd, check=True)


movie_args = '--movie-time 10^(1/8)'.split()


def run_replicas(name, max_iter=max_iter_default, min_T=0.001, seed=None, max_independent_samples=None, extraname='', extraflags=''):
    save = f'z+{name}+seed-{seed}+{extraname}'
    samples = []
    seed_str = []
    if seed is not None:
        seed_str = ['--seed', str(seed)]
    if max_independent_samples is not None:
        samples = ['--max-independent-samples', str(max_independent_samples)]
    rq(name=save,
       cmd=['../target/release/replicas']+systems[name]+movie_args
        + f'--save-time 0.5 --save-as {save}.cbor'.split()
        + extraflags.split()
        + f'--max-iter {max_iter} --min-T {min_T}'.split()
        + samples
        + seed_str,
       cpus='all')

def geometric_spacing(min, max, number):
    numbers = []
    r = (max/min)**(1/(number-1))
    for i in range(number):
        numbers.append(min*r**i)
    return numbers

def run_tempering(name, max_iter=max_iter_default, min_T=0.001, max_T=1, num_T=20, can_steps=1, seed=None, max_independent_samples=None, extraname='', extraflags=''):
    T = geometric_spacing(min_T,max_T,num_T)
    T_string = ''
    seed_str = []
    if seed is not None:
        seed_str = ['--seed', str(seed)]
    for t in T:
        T_string += '--T '+str(t)+' '
    save = f'tem+{name}{extraname}+seed-{seed}'
    rq(name=save,
       cmd=['../target/release/tempering']+systems[name]+movie_args
        + f'--save-time 0.5 --save-as {save}.cbor'.split()
        + extraflags.split()
        + f'--max-iter {max_iter} {T_string}'.split()
        + f'--canonical-steps {can_steps}'.split()
        + seed_str,
       cpus='all')


def histogram(name, de, translation_scale, seed_str):
    return f'../target/release/histogram --save-time 0.5 --energy-bin {de} --translation-scale {translation_scale} {seed_str}'.split()+movie_args+systems[name]


def run_sad(name, de, max_iter=max_iter_default, min_T=0.001, max_E=None, translation_scale=0.05, seed=None, extraname=''):
    de = str(de)
    if seed is None:
        seed_str = ''
        save = f'sad-{name}-{de}+{translation_scale}'
    else:
        seed_str = f'--seed {seed}'
        save = f'sad+{name}+seed-{seed}+de-{de}+step-{translation_scale}'
    max_E_args = []
    if max_E is not None:
        max_E_args = f'--max-allowed-energy {max_E}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale, seed_str=seed_str)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --sad-min-T {min_T}'.split()
        + max_E_args,
       cpus=1)

#TODO: Make it accept a seed
def run_wl(name, de, min_E, max_E, min_gamma=None, max_iter=max_iter_default, translation_scale=0.05, seed=None):
    de = str(de)
    if seed is None:
        seed_str = ''
        save = f'wl+{name}-{de}+{translation_scale}'
    else:
        seed_str = f'--seed {seed}'
        save = f'wl+{name}+seed-{seed}+de-{de}+step-{translation_scale}'
    min_gamma_args = []
    if min_gamma is not None:
        min_gamma_args = f'--wl-min-gamma {min_gamma}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale, seed_str=seed_str)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split()
        + min_gamma_args,
       cpus=1)


def run_inv_t_wl(name, de, min_E, max_E, max_iter=max_iter_default, translation_scale=0.05, seed=None):
    de = str(de)
    if seed is None:
        seed_str = ''
        save = f'itwl+{name}-{de}+{translation_scale}'
    else:
        seed_str = f'--seed {seed}'
        save = f'itwl+{name}+seed-{seed}+de-{de}+step-{translation_scale}'
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale, seed_str=seed_str)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --inv-t-wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split(),
       cpus=1)


min_T = 0.0001

E1 = -133.58642  # minimum energy (Mackay) for an LJ31 cluster
E2 = -133.29382  # first local minimum (anti-Mackay) for an LJ31 cluster
E_transition = -131  # approximate energy of the transition state between the two
T_transition = 0.025  # approximate temperature for transition between the two

T_trans_1_r2 = system.systems['T-trans-1']['R_small']
T_trans_1_h2 = system.systems['T-trans-1']['h_small']
T_trans_1_n = system.systems['T-trans-1']['n']

hard_r2 = system.systems['hard']['R_small']
hard_h2 = system.systems['hard']['h_small']
hard_n = system.systems['hard']['n']

tiny_r2 = system.systems['tiny']['R_small']
tiny_h2 = system.systems['tiny']['h_small']
tiny_n = system.systems['tiny']['n']

systems = {
    'lj31-like': '--two-wells-N 90 --two-wells-h2-to-h1 1.005 --two-wells-barrier-over-h1 0.03 --two-wells-r2 0.75'.split(),

    #For the thesis
    'T-trans-1+barrier-1': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 1.0 --two-wells-r2 {T_trans_1_r2}'.split(),
    'T-trans-1+barrier-4e-1': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 0.4 --two-wells-r2 {T_trans_1_r2}'.split(),
    'T-trans-1+barrier-2e-1': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 0.2 --two-wells-r2 {T_trans_1_r2}'.split(),
    'T-trans-1+barrier-1e-1': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 0.1 --two-wells-r2 {T_trans_1_r2}'.split(),
    'T-trans-1+barrier-0': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 0 --two-wells-r2 {T_trans_1_r2}'.split(),

    'T-trans-1+barrier-1-small': f'--two-wells-N {3} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 1.0 --two-wells-r2 {T_trans_1_r2}'.split(),

    'easier': '--two-wells-N 9 --two-wells-h2-to-h1 1.1347 --two-wells-barrier-over-h1 0.5 --two-wells-r2 0.5'.split(),
    'easier-all-barrier': '--two-wells-N 9 --two-wells-h2-to-h1 1.1347 --two-wells-barrier-over-h1 1 --two-wells-r2 0.5'.split(),
    'easier-no-barrier': '--two-wells-N 9 --two-wells-h2-to-h1 1.1347 --two-wells-barrier-over-h1 0 --two-wells-r2 0.5'.split(),

    'easiest': '--two-wells-N 9 --two-wells-h2-to-h1 1 --two-wells-barrier-over-h1 0.5 --two-wells-r2 0.5'.split(),

    'tiny': f'--two-wells-N {tiny_n} --two-wells-h2-to-h1 {tiny_h2} --two-wells-barrier-over-h1 0 --two-wells-r2 {tiny_r2}'.split(),

    'hard': f'--two-wells-N {hard_n} --two-wells-h2-to-h1 {hard_h2} --two-wells-barrier-over-h1 0.1 --two-wells-r2 {hard_r2}'.split(),
    'hard-half-barrier': f'--two-wells-N {hard_n} --two-wells-h2-to-h1 {hard_h2} --two-wells-barrier-over-h1 0.05 --two-wells-r2 {hard_r2}'.split(),
    'hard-fifth-barrier': f'--two-wells-N {hard_n} --two-wells-h2-to-h1 {hard_h2} --two-wells-barrier-over-h1 0.02 --two-wells-r2 {hard_r2}'.split(),
    'hard-no-barrier': f'--two-wells-N {hard_n} --two-wells-h2-to-h1 {hard_h2} --two-wells-barrier-over-h1 0 --two-wells-r2 {hard_r2}'.split(),
}
seeds = [1,12,123,1234,12345,123456,1234567,12345678]
#run_replicas(name='tiny', min_T=system.systems['tiny']['min_T'], max_iter=1e13, max_independent_samples=3e7,
#             extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

## Uncomment when you want to rerun stuff ##

# for seed in seeds:
#     for s in ['T-trans-1+barrier-2e-1']:
#         for de in [1e-5,1e-4]:#, 0.0001]:
#             for translation_scale in [1e-2,1e-3,1e-4]:#[1e-4, 1e-3, 0.01]:
#                 if de == 1e-4 and translation_scale != 1e-3:
#                     pass
#                 else:

# #                     run_sad(name=s, min_T=system.systems['T-trans-1']['min_T'], max_iter=1e13, translation_scale=translation_scale, de=de, seed=seed)
                
#                 ## UNCOMMENT THESE WHEN YOU KNOW WHICH TRANSLATION SCALE TO USE ##


#                     run_inv_t_wl(name=s, min_E=system.systems['T-trans-1']['min_E'], max_E=de/2, max_iter=1e12,
#                             translation_scale=translation_scale, de=de, seed=seed)

for seed in seeds:
    for s in ['T-trans-1+barrier-0', 'T-trans-1+barrier-1e-1', 'T-trans-1+barrier-2e-1']:
        for de in [1e-5]:#, 0.0001]:
            for translation_scale in [1e-4, 1e-3, 0.01]:
                if de == 1e-4 and translation_scale != 1e-3:
                    pass
                else:
                    pass
                    #run_sad(name=s, min_T=system.systems['T-trans-1']['min_T'], max_iter=1e12, translation_scale=translation_scale, de=de, seed=seed)
                
                ## UNCOMMENT THESE WHEN YOU KNOW WHICH TRANSLATION SCALE TO USE ##

                    #run_inv_t_wl(name=s, min_E=system.systems['T-trans-1']['min_E'], max_E=de/2, max_iter=1e12,
                    #        translation_scale=translation_scale, de=de, seed=seed)



#run_tempering('T-trans-1+barrier-0', max_iter=1e12, num_T=10, can_steps=10, extraname='temps-10-')

# for seed in seeds:
#     for s in ['T-trans-1+barrier-0', 'T-trans-1+barrier-1e-1', 'T-trans-1+barrier-2e-1']:
#         run_tempering(s, max_iter=1e13, num_T=50, can_steps=10, seed=seed)

seed = 1
for s in ['T-trans-1+barrier-1', 'T-trans-1+barrier-1-small']:
    run_replicas(s, min_T=system.systems['T-trans-1']['min_T'], max_iter=1e13, max_independent_samples=3e7,
             extraflags=' --independent-systems-before-new-bin 16', seed=seed)
