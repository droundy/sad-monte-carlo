_markers = {
    '0.01+0.01': 'D', '0.01+0.001': '^', '0.001+0.01': 'o', '0.001+0.001': 'x',
    '0.001+0.1': '^', '0.01+0.1': 'x',
    '0.0001+0.01': '^', '0.0001+0.001': '<', '0.0001+0.0001': '>',
    '05+0.01': 'x', '05+0.001': 'D', '05+0.0001': 'o',
}

_colors = {'z': 'k', 'wl': 'b', 'itwl': 'g', 'sad': 'tab:orange'}
dashes = {0: 'solid', 1: 'dashed'}
_linestyles = {
    'z': '-',
    'wl': 'dashed',
    'itwl': 'dashdot',
    'sad': (0, (3, 1, 1, 1, 1, 1))
}


def marker(base):
    if base is None:
        return None
    precision = base[base.rfind('-') + 1:]
    if precision in _markers:
        return _markers[precision]
    else:
        return None


def color(base):
    if base is None:
        return 'tab:cyan'
    method = base[:base.find('-')]
    if method in _colors:
        return _colors[method]
    else:
        return None


def linestyle(base):
    if base is None:
        return ':'
    method = base[:base.find('-')]
    if method in _linestyles:
        return _linestyles[method]
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
            dE = r'1\times 10^5'

        ds = base[base.find('+')+1:]
        if method_name == 'itwl':
            name = '$1/t$-WL'
        else:
            name = method_name.upper()
        # for proposal, we will omit any parameters that are identical for all plots
        return rf'{name}, $\Delta E={dE}$'
        return name # for if we use same delta E and translation scale for each
        return rf'{name}, $\Delta E={dE}$, $\Delta s={ds}$'
    except:
        return base