def readsystem(data):
    system = {
        'kind': 'unknown',
    }
    data_system = None
    if 'system' in data:
        data_system = data['system']
    else:
        data_system = data['system_with_lowest_max_energy'][0]
    if 'Wca' in data_system:
        system['kind'] = 'WCA'
        a = data_system['Wca']['cell']['box_diagonal']
        system['volume'] = a['x']*a['y']*a['z']
        system['N'] = len(data_system['Wca']['cell']['positions'])
        system['density'] = system['N']/system['volume']
    elif 'Fake' in data_system:
        fake = data_system['Fake']
        if 'Gaussian' in fake['function']:
            system['kind'] = 'Gaussian'
            system['sigma'] = fake['function']['Gaussian']['sigma']
        elif 'Pieces' in fake['function']:
            system['kind'] = 'Pieces'
            system['a'] = fake['function']['Pieces']['a']
            system['b'] = fake['function']['Pieces']['b']
            system['e1'] = fake['function']['Pieces']['e1']
            system['e2'] = fake['function']['Pieces']['e2']
    elif 'FakeErfinv' in data_system:
        system['kind'] = 'Erfinv'
        system['mean_energy'] = data_system['FakeErfinv']['parameters']['mean_energy']
        system['N'] = len(data_system['FakeErfinv']['position'])
    return system
