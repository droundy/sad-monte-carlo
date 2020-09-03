def readsystem(data):
    system = {}
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
    return system