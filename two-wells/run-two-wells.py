#!/usr/bin/python3

import numpy as np
from subprocess import run
import system

run(['cargo', 'build', '--release', '--bin',
     'replicas', '--bin', 'histogram'], check=True)

max_iter_default = 1e13


def rq(name, cmd, cpus):
    run(f'rq run -c {cpus} --max-output=30 -R -J'.split() +
        [name, '--']+cmd, check=True)


movie_args = '--movie-time 10^(1/8)'.split()


def run_replicas(name, max_iter=max_iter_default, min_T=0.001, max_independent_samples=None, extraname='', extraflags=''):
    save = f'z-{extraname}{name}'
    samples = []
    if max_independent_samples is not None:
        samples = ['--max-independent-samples', str(max_independent_samples)]
    rq(name=save,
       cmd=['../target/release/replicas']+systems[name]+movie_args
        + f'--save-time 0.5 --save-as {save}.cbor'.split()
        + extraflags.split()
        + f'--max-iter {max_iter} --min-T {min_T}'.split()
        + samples,
       cpus='all')


def histogram(name, de, translation_scale, seed_str):
    return f'../target/release/histogram --save-time 0.5 --energy-bin {de} --translation-scale {translation_scale} {seed_str}'.split()+movie_args+systems[name]


def run_sad(name, de, max_iter=max_iter_default, min_T=0.001, max_E=None, translation_scale=0.05, seed=None, extraname=''):
    de = str(de)
    if seed is not None:
        seed_str = ''
        save = f'sad-{name}-{de}+{translation_scale}'
    else:
        seed_str = f'--seed {seed}'
        save = f'sad-{name}-{seed}-{de}+{translation_scale}'
    max_E_args = []
    if max_E is not None:
        max_E_args = f'--max-allowed-energy {max_E}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale, seed = seed_str)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --sad-min-T {min_T}'.split()
        + max_E_args,
       cpus=1)

#TODO: Make it accept a seed
def run_wl(name, de, min_E, max_E, min_gamma=None, max_iter=max_iter_default, translation_scale=0.05):
    de = str(de)
    save = f'wl-{name}-{de}+{translation_scale}'
    min_gamma_args = []
    if min_gamma is not None:
        min_gamma_args = f'--wl-min-gamma {min_gamma}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split()
        + min_gamma_args,
       cpus=1)


def run_inv_t_wl(name, de, min_E, max_E, max_iter=max_iter_default, translation_scale=0.05):
    de = str(de)
    save = f'itwl-{name}-{de}+{translation_scale}'
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale)
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
    'T-trans-1+1e-1': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 0.1 --two-wells-r2 {T_trans_1_r2}'.split(),
    'T-trans-1+0': f'--two-wells-N {T_trans_1_n} --two-wells-h2-to-h1 {T_trans_1_h2} --two-wells-barrier-over-h1 0 --two-wells-r2 {T_trans_1_r2}'.split(),

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

for seed in seeds:
    for s in ['T-trans-1+0', 'T-trans-1+1e-1']:
        for de in [0.00001, 0.0001]:
            for translation_scale in [1e-4, 1e-3, 0.01]:
                    run_sad(name=s, min_T=system.systems['T_trans_1']['min_T'], max_iter=1e12, translation_scale=translation_scale, de=de, seed=seed)
                
                ## UNCOMMENT THESE WHEN YOU KNOW WHICH TRANSLATION SCALE TO USE ##

                #run_wl(name=s, min_E=system.systems[s]['min_E'], max_E=de/2, max_iter=1e12,
                #            translation_scale=translation_scale, de=de, min_gamma=1e-9)
                #run_inv_t_wl(name=s, min_E=system.systems[s]['min_E'], max_E=de/2, max_iter=1e12,
                #            translation_scale=translation_scale, de=de)


# hard_min_T = system.systems['hard']['min_T']
# hard_min_E = system.systems['hard']['min_E']
# run_replicas(name='hard-half-barrier', min_T=hard_min_T, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
# run_replicas(name='hard-fifth-barrier', min_T=hard_min_T, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

# for s in ['hard-half-barrier','hard-no-barrier', 'hard']:#['hard', 'hard-no-barrier', 'hard-fifth-barrier']:
#     for de in [0.001, 0.01]:
#         for translation_scale in [0.01, 0.1]:
#                 run_wl(name=s, min_E=hard_min_E, max_E=de/2, max_iter=1e12,
#                             translation_scale=translation_scale, de=de, min_gamma=1e-9)
#                 run_inv_t_wl(name=s, min_E=hard_min_E, max_E=de/2, max_iter=1e12,
#                             translation_scale=translation_scale, de=de)
#                 run_sad(name=s, min_T=hard_min_T, max_iter=1e12, translation_scale=translation_scale, de=de)

# run_replicas(name='hard-no-barrier', min_T=hard_min_T, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
# run_replicas(name='hard', min_T=hard_min_T, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

# run_replicas(name='easier-all-barrier', min_T=0.005, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
# run_replicas(name='easier', min_T=0.005, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
# run_replicas(name='easier-no-barrier', min_T=0.005, max_iter=1e13, max_independent_samples=10000,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

# for de in [0.001, 0.01, 0.1]:
#     for translation_scale in [0.001, 0.01]:
#         for system in ['easier-no-barrier', 'easier-all-barrier', 'easier']:
#             run_wl(name=system, min_E=-1.13, max_E=0.0005, max_iter=1e12,
#                         translation_scale=translation_scale, de=de, min_gamma=1e-9)
#             run_inv_t_wl(name=system, min_E=-1.13, max_E=0.0005, max_iter=1e12,
#                         translation_scale=translation_scale, de=de)
#             run_sad(name=system, min_T=0.005, max_iter=1e12, translation_scale=translation_scale, de=de)

# run_inv_t_wl(name='easiest', min_E=-1.13, max_E=0.0005, max_iter=1e12,
#              translation_scale=0.05, de=0.001)
# run_replicas(name='easiest', min_T=0.005, max_iter=1e13, max_independent_samples=100000)

# run_replicas(name='hard', min_T=0.001, max_iter=1e13, max_independent_samples=100)
# run_replicas(name='hard', min_T=0.001, max_iter=1e13, max_independent_samples=100,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

# run_replicas(name='lj31-like', min_T=min_T, max_iter=1e13, max_independent_samples=100,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
# run_replicas(name='lj31-like', min_T=min_T,
#              max_iter=1e13, max_independent_samples=100)
