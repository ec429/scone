#!/usr/bin/python

import math
import cfg

with open('rssk.cfg') as f:
    kop = cfg.parse(f)

kc = kop['@Kopernicus[Kerbol?System]:FOR[RealSolarSystem]'][0]

bodies = kc['Body']

sun = [b for b in bodies if b['name'] == 'Sun'][0]
sungp = float(sun['Properties'][0]['gravParameter'])
print "double sungp = %g;"%(sungp,)

planets = [b for b in bodies if b.get('Orbit', [{}])[0].get('referenceBody') == 'Sun']

print """struct planet {
\tconst char *name;
\tdouble sma;
\tdouble ecc;
\tdouble inc;
\tdouble lan;
\tdouble ape;
\tdouble mar;
} planets[] = {"""
for p in planets:
    if 'cbNameLater' in p:
        p['name'] = p.pop('cbNameLater')
    o = p['Orbit'][0]
    data = {'name': p['name'],
            'sma': float(o['semiMajorAxis']),
            'ecc': float(o['eccentricity']),
            'inc': float(o['inclination']),
            'lan': float(o['longitudeOfAscendingNode']),
            'ape': float(o['argumentOfPeriapsis']),
            'mar': float(o['meanAnomalyAtEpochD']) * math.pi / 180.0,
            }
    print "\t{"
    for k,v in data.items():
        if isinstance(v, str):
            v = '"' + v + '"'
        print "\t\t.%s = %s,"%(k, v)
    print "\t},"
print "};"
