#!/usr/bin/python3

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

def run_sad(name, de, max_iter=max_iter_default, min_T=0.001, max_E=None, translation_scale=0.05):
    de = str(de)
    save = 'sad-'+name+'-'+de
    max_E_args = []
    if max_E is not None:
        max_E_args = f'--max-allowed-energy {max_E}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --sad-min-T {min_T}'.split()
        + max_E_args,
       cpus='1')


def run_wl(name, de, min_E, max_E, min_gamma=None, max_iter=max_iter_default, translation_scale=0.05):
    de = str(de)
    save = 'wl-'+name+'-'+de
    min_gamma_args = []
    if min_gamma is not None:
        min_gamma_args = f'--wl-min-gamma {min_gamma}'.split()
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split()
        + min_gamma_args,
       cpus='1')


def run_inv_t_wl(name, de, min_E, max_E, max_iter=max_iter_default, translation_scale=0.05):
    de = str(de)
    save = 'itwl-'+name+'-'+de
    rq(name=save,
       cmd=histogram(name, de, translation_scale=translation_scale)
        + f'--save-as {save}.cbor'.split()
        + f'--max-iter {max_iter} --inv-t-wl --min-allowed-energy {min_E} --max-allowed-energy {max_E}'.split(),
       cpus='1')


systems = {
    'lj31':  '--lj-N 31 --lj-radius 2.5 --independent-systems-before-new-bin 16'.split(),
    'big-lj31':  '--lj-N 31 --lj-radius 5 --independent-systems-before-new-bin 16'.split(),
    'huge-lj31':  '--lj-N 31 --lj-radius 15 --independent-systems-before-new-bin 16'.split(),
}

run_replicas(name='lj31', min_T=0.001, max_iter=1e14, max_independent_samples=1000)
run_replicas(name='big-lj31', min_T=0.001, max_iter=1e14, max_independent_samples=1000)
run_replicas(name='huge-lj31', min_T=0.001, max_iter=1e14, max_independent_samples=1000)
exit(1)
run_replicas(name='huge-lj31', min_T=0.001, max_iter=1e14)
run_replicas(name='lj31', min_T=0.001, extraname='0.001-', max_iter=1e14)
run_replicas(name='lj31', min_T=0.005, max_iter=1e14)
# run_replicas(name='lj31', min_T=0.005, extraname='mean-', extraflags='--mean-for-median')
run_replicas(name='biglj31', min_T=0.005, max_iter=1e12)

for de in [0.01]: # , 1/2**10]: 0.1
    run_sad('lj31', de=de, min_T=0.005, max_E=0)
    # run_wl('lj31', de=de, min_gamma=1e-10, min_E=-133.53, max_E=-110)
    run_wl('lj31', de=de, min_gamma=1e-10, min_E=-133.53, max_E=0)
    run_inv_t_wl('lj31', de=de, min_E=-133.53, max_E=0)
