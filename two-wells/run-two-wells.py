#!/usr/bin/python3

import numpy as np
from subprocess import run

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


def histogram(name, de, translation_scale):
    return f'../target/release/histogram --save-time 0.5 --energy-bin {de} --translation-scale {translation_scale}'.split()+movie_args+systems[name]


def run_sad(name, de, max_iter=max_iter_default, min_T=0.001, max_E=None, translation_scale=0.05, extraname=''):
    de = str(de)
    save = f'sad-{name}-{de}+{translation_scale}'
    max_E_args = []
    if max_E is not None:
        max_E_args = f'--max-allowed-energy {max_E}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --sad-min-T {min_T}'.split()
        + max_E_args,
       cpus=1)


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

systems = {
    'lj31-like': '--two-wells-N 90 --two-wells-h2-to-h1 1.005 --two-wells-barrier-over-h1 0.03 --two-wells-r2 0.75'.split(),
    'easy': '--two-wells-N 30 --two-wells-h2-to-h1 1.1 --two-wells-barrier-over-h1 0.5 --two-wells-r2 0.5'.split(),
    'easier': '--two-wells-N 9 --two-wells-h2-to-h1 1.1347 --two-wells-barrier-over-h1 0.5 --two-wells-r2 0.5'.split(),
    'easier-all-barrier': '--two-wells-N 9 --two-wells-h2-to-h1 1.1347 --two-wells-barrier-over-h1 1 --two-wells-r2 0.5'.split(),
    'easier-no-barrier': '--two-wells-N 9 --two-wells-h2-to-h1 1.1347 --two-wells-barrier-over-h1 0 --two-wells-r2 0.5'.split(),
    'easiest': '--two-wells-N 9 --two-wells-h2-to-h1 1 --two-wells-barrier-over-h1 0.5 --two-wells-r2 0.5'.split(),
}

run_replicas(name='easier-all-barrier', min_T=0.005, max_iter=1e13, max_independent_samples=10000,
             extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
run_replicas(name='easier', min_T=0.005, max_iter=1e13, max_independent_samples=10000,
             extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
run_replicas(name='easier-no-barrier', min_T=0.005, max_iter=1e13, max_independent_samples=10000,
             extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

for translation_scale in [0.001, 0.01]:
    for system in ['easier-no-barrier', 'easier-all-barrier', 'easier']:
        run_wl(name=system, min_E=-1.13, max_E=0.0005, max_iter=1e12,
                    translation_scale=translation_scale, de=0.001, min_gamma=1e-9)
        run_inv_t_wl(name=system, min_E=-1.13, max_E=0.0005, max_iter=1e12,
                    translation_scale=translation_scale, de=0.001)
        run_sad(name=system, min_T=0.005, max_iter=1e12, translation_scale=translation_scale, de=0.001)

# run_inv_t_wl(name='easiest', min_E=-1.13, max_E=0.0005, max_iter=1e12,
#              translation_scale=0.05, de=0.001)
# run_replicas(name='easiest', min_T=0.005, max_iter=1e13, max_independent_samples=100000)

# run_replicas(name='easy', min_T=0.001, max_iter=1e13, max_independent_samples=100)
# run_replicas(name='easy', min_T=0.001, max_iter=1e13, max_independent_samples=100,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')

# run_replicas(name='lj31-like', min_T=min_T, max_iter=1e13, max_independent_samples=100,
#              extraflags=' --independent-systems-before-new-bin 16', extraname='i16-')
# run_replicas(name='lj31-like', min_T=min_T,
#              max_iter=1e13, max_independent_samples=100)
