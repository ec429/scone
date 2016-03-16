#!/usr/bin/python

def parse(f):
    top = {}
    current = top
    stack = []
    last = None

    for line in f:
        line = line.strip()
        if line.startswith('//'): continue
        if '=' in line:
            k, _, v = line.partition('=')
            current[k.strip()] = v.strip()
            continue
        if line == '{':
            new = {}
            current.setdefault(last, []).append(new)
            stack.append(current)
            current = new
            continue
        if line == '}':
            current = stack.pop()
            last = None
            continue
        last = line
    return(top)

if __name__ == '__main__':
    import sys, pprint
    d = parse(sys.stdin)
    pprint.pprint(d)
