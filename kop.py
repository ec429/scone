#!/usr/bin/python

from math import radians
import cfg

with open('rssk.cfg') as f:
    kop = cfg.parse(f)

kc = kop['@Kopernicus[Kerbol?System]:FOR[RealSolarSystem]'][0]

print "double epoch = %g;"%(float(kc['Epoch']),)

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
            'inc': radians(float(o['inclination'])),
            'lan': radians(float(o['longitudeOfAscendingNode'])),
            'ape': radians(float(o['argumentOfPeriapsis'])),
            'mar': radians(float(o['meanAnomalyAtEpochD'])),
            }
    print "\t{"
    for k,v in data.items():
        if isinstance(v, str):
            v = '"' + v + '"'
        print "\t\t.%s = %s,"%(k, v)
    print "\t},"
print "};"
