_markers= {'0.01+0.01':'D','0.01+0.001':'^','0.001+0.01':'o','0.001+0.001':'x',
            '0.001+0.1':'^','0.01+0.1':'x',
            '05+0.001':'x',
            '05+0.0001':'x',
            '05+0.00001':'x',
            '0.0001+0.01':'^','1e-05+0.1':'x','0.0001+0.1':'o','1e-05+0.01':'D','05+0.1':'x','05+0.01':'D'}

_colors = {'z':'k','wl':'b','itwl':'g','sad':'tab:orange'}
dashes = {0: 'solid', 1:'dashed'}
_linestyles = {'z':'solid','wl':'dashed','itwl':'dashdot','sad':(0, (3, 1, 1, 1, 1, 1))}

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
        return None
    method = base[:base.find('-')]
    if method in _colors:
        return _colors[method]
    else:
        return None

def linestyle(base):
    if base is None:
        return None
    method = base[:base.find('-')]
    if method in _linestyles:
        return _linestyles[method]
    else:
        return None
