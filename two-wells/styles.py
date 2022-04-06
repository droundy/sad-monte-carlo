from find_phase_transition import S


_markers = {
    'de-0.0001+step-0.001': 'D', 'de-0.0001+step-0.01': '^', 'de-0.0001+step-0.0001': 'o', 'de-0.0001+0.001': 'x',
    'de-0.0001+0.1': '^', 'de-1e-05+step-0.001': 'x',
    'sad': 's', 'itwl': '<', 'z': 'o', 'tem': 'x',
    'de-1e-05+step-0.001': '>',
    'de-1e-05+step-0.1': 'x', '05+0.001': 'D', '05+0.0001': 'o',
}

#_colors = {'z': 'k', 'wl': 'b', 'itwl': 'g', 'sad': 'tab:orange'}
_colors = {'z': 'k', 'wl': 'b', 'barrier-0': 'g', 'barrier-1e-1': 'tab:orange', 'barrier-2e-1': 'tab:cyan'}
dashes = {0: 'solid', 1: 'dashed'}
_linestyles = {
    'barrier-0': '-',
    'barrier-1e-1': 'dashed',
    'barrier-2e-1': '-.',
    'itwl': 'dashdot',
    'sad': (0, (3, 1, 1, 1, 1, 1))
}


# def marker(base):
#     if base is None:
#         return None
#     # return None # Always use no marker if dE etc. is same for all simulations
#     splits = base.split('+')
#     precision = splits[-2] + '+' + splits[-1]
#     if precision in _markers:
#         return _markers[precision]
#     else:
#         return None

def marker(base):
    if base is None:
        return None
    if 'itwl' in base:
        return _markers['itwl']
    if 'sad' in base:
        return _markers['sad']
    if 'z' in base:
        return _markers['z']
    if 'tem' in base:
        return _markers['tem']


def color(base):
    if base is None:
        return 'b'
    barrier = 'barrier-'+get_barrier(base)
    if barrier in _linestyles:
        return _colors[barrier]
    else:
        return None
    # if base is None:
    #     return ':'
    # barr_ind = base.find('barrier')
    # barrier = base[barr_ind:barr_ind+base[barr_ind:].find('+')]
    # if barrier in _colors:
    #     return _colors[barrier]
    # else:
    #     return None
    # if base is None:
    #     return 'tab:cyan'
    # method = base[:base.find('+')]
    # if method in _colors:
    #     return _colors[method]
    # else:
    #     return None

def get_barrier(base):
    if base is None:
        return 'exact'
    b = base.split('barrier-')[-1]
    return b.split('+')[0]


def linestyle(base):
    if base is None:
        return ':'
    barrier = 'barrier-'+get_barrier(base)
    if barrier in _linestyles:
        return _linestyles[barrier]
    else:
        return None

def pretty_label(base):
    if base is None:
        return 'exact'
    if 'z-' in base:
        return 'ZMC'
    
    try:
        method_name = base[:base.find('-')]
        #find index of first digit
        for i, e in enumerate(base):
            if e.isdigit():
                dE = base[i:base.find('+')]
                break
        if dE == '1e-05':
            dE = r'10^{-5}'

        ds = base[base.find('+')+1:]
        if method_name == 'itwl':
            name = '$1/t$-WL'
        else:
            name = method_name.upper()
        # for proposal, we will omit any parameters that are identical for all plots
        # return rf'{name}, $\Delta E={dE}$'
        return name # for if we use same delta E and translation scale for each
        return rf'{name}, $\Delta E={dE}$, $\Delta s={ds}$'
    except:
        return base