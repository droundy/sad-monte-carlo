def readsystem(data):
    system = {
        'kind': 'unknown',
    }
    if 'Wca' in data['system']:
        system['kind'] = 'WCA'
        a = data['system']['Wca']['cell']['box_diagonal']
        system['volume'] = a['x']*a['y']*a['z']
        system['N'] = len(data['system']['Wca']['cell']['positions'])
        system['density'] = system['N']/system['volume']
    elif 'Fake' in data['system']:
        fake = data['system']['Fake']
        if 'Gaussian' in fake['function']:
            system['kind'] = 'Gaussian'
            system['sigma'] = fake['function']['Gaussian']['sigma']
        elif 'Pieces' in fake['function']:
            system['kind'] = 'Pieces'
            system['a'] = fake['function']['Pieces']['a']
            system['b'] = fake['function']['Pieces']['b']
            system['e1'] = fake['function']['Pieces']['e1']
            system['e2'] = fake['function']['Pieces']['e2']
    elif 'FakeErfinv' in data['system']:
        system['kind'] = 'Erfinv'
        system['mean_energy'] = data['system']['FakeErfinv']['parameters']['mean_energy']
        system['N'] = len(data['system']['FakeErfinv']['position'])
    return system
