# Genes we want to plot
genes_of_interest = {
    'HNEAP_RS07320': dict(name='Crp/Fnr', kind='regulation'),
    'HNEAP_RS01030': dict(name='DAB2B', kind='transport'),
    'HNEAP_RS01035': dict(name='DAB2A', kind='transport'), 
    'HNEAP_RS04585': dict(name='DAB1B', kind='transport'),
    'HNEAP_RS04595': dict(name='DAB1A', kind='transport'),
    'HNEAP_RS01040': dict(name='LysR DAB2', kind='regulation'),
    'HNEAP_RS04645': dict(name='csosCA', kind='CA'),
    'HNEAP_RS04565': dict(name='csos1D', kind='carboxysome'),
    'HNEAP_RS04655': dict(name='cbbS', kind='carboxysomal rubisco'),
    'HNEAP_RS04660': dict(name='cbbL', kind='carboxysomal rubisco'),
    'HNEAP_RS04615': dict(name='acRAF', kind='chaperone'),
    'HNEAP_RS04620': dict(name='csos1B', kind='carboxysome'),
    'HNEAP_RS04625': dict(name='csos1A', kind='carboxysome'),
    'HNEAP_RS05490': dict(name='LysR', kind='regulation'),
    'HNEAP_RS04640': dict(name='csos4A', kind='carboxysome'),
    'HNEAP_RS04635': dict(name='csos4B', kind='carboxysome'),
    'HNEAP_RS04650': dict(name='csos2', kind='carboxysome'),
    'HNEAP_RS04630': dict(name='csos1C', kind='carboxysome'), 
    # Non-carboxysomal FII rubisco was named cbbM in Baker et al. JBac 1998.
    'HNEAP_RS05505': dict(name='cbbM', kind='non-carboxysomal rubisco'),
    'HNEAP_RS04600': dict(name='cbbO', kind='chaperone'),
    'HNEAP_RS04575': dict(name='cbbQ', kind='chaperone'),
    # mcdAB from MacCready 2021
    'HNEAP_RS04610': dict(name='mcdA', kind='regulation'),
    'HNEAP_RS12660': dict(name='mcdB', kind='regulation'),
}

categories2plot = {
    'transport': True,
    'carboxysome': True,
    'CA': True,
    'carboxysomal rubisco': True,
    'regulation': False,
    'chaperone': False,
    'non-carboxysomal rubisco': False
}

genes_of_interest_filtered = dict((k, v) for (k,v) in genes_of_interest.items()
                                  if categories2plot[v['kind']])