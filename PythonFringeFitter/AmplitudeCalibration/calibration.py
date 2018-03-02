importfitsidi(fitsidifile=['n14c3_1_1.IDI1', 'n14c3_1_1.IDI2'],vis='n14c3.ms',constobsid=True,scanreindexgap_s=15)

execfile('tsys_spectrum.py')
gencal(vis='n14c3.ms', caltable='n14c3.tsys', caltype='tsys')

execfile('fixms.py')
applycal(vis='n14c3.ms', gaintable=['n14c3.tsys', 'n14c3.gc'])
