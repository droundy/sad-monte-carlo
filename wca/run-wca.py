#!/usr/bin/python3

from subprocess import run
import numpy as np

run(['cargo', 'build', '--release', '--bin',
     'replicas', '--bin', 'binning'], check=True)

max_iter_default = 1e11

def rq(name, cmd, cpus):
    run(f'rq run -c {cpus} --max-output=30 -R -J'.split() +
        [name, '--']+cmd, check=True)


movie_args = '--movie-time 10^(1/4)'.split()

def run_replicas(name, max_iter=max_iter_default, min_T=0.001):
    save = 'r-'+name
    rq(name=save,
       cmd=['../target/release/replicas']+systems[name]+movie_args
        + f' --save-time 0.5 --save-as {save}.yaml'.split()
        + f'--max-iter {max_iter} --min-T {min_T}'.split(),
       cpus='all')

def binning_histogram(name, de, translation_scale):
    return f'../target/release/binning --save-time 0.5 --histogram-bin {de} --translation-scale {translation_scale}'.split()+movie_args+systems[name]

def run_sad(name, de, max_iter=max_iter_default, min_T=0.001, translation_scale=0.05):
    de = str(de)
    save = 'sad-'+name+'-'+de
    rq(name=save,
       cmd=binning_histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.yaml'.split()
        + f'--max-iter {max_iter} --sad-min-T {min_T}'.split(),
       cpus='1')


def run_wl(name, de, min_E, max_E, max_iter=max_iter_default, translation_scale=0.05):
    de = str(de)
    save = 'wl-'+name+'-'+de
    rq(name=save,
       cmd=binning_histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.yaml'.split()
        + f'--max-iter {max_iter} --wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split(),
       cpus='1')


def run_inv_t_wl(name, de, min_E, max_E, max_iter=max_iter_default, translation_scale=0.05):
    de = str(de)
    save = 'itwl-'+name+'-'+de
    rq(name=save,
       cmd=binning_histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.yaml'.split()
        + f'--max-iter {max_iter} --inv-t-wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split(),
       cpus='1')

systems = {
    'pieces': '--fake-pieces-a 0.1 --fake-pieces-b 0.2 --fake-pieces-e1 1.0 --fake-pieces-e2 0.5'.split(),
    'erfinv': '--fake-erfinv-mean-energy 0 --fake-erfinv-N 3'.split(),
    'linear': '--fake-linear'.split(),
    'quadratic': '--fake-quadratic-dimensions 3'.split(),
}

densities = np.arange(0.75, 1.4, 0.125)
min_T = 0.5

for d in densities:
    name = f'wca-32-{d}'
    systems[name] = f'--wca-reduced-density {d} --wca-N 32'
    run_replicas(name=name, min_T = min_T)
