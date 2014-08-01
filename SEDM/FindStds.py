import pyfits as pf
import NPK.Standards as SS
import sys

stds = SS.Standards.keys()
stdset = map(set, stds)

out = "["
if __name__ == '__main__':
    
    files = sys.argv[1:]

    for file in files:
        FF = pf.open(file)
        name = FF[0].header['NAME'].lower().lstrip('std-')
        name = name.replace('+','')
        
        std = False
        for i in xrange(len(stds)):
            std_name = stdset[i]
            if set(name).issubset(std_name):
                std = stds[i]

        if std:
            out += '{"infiles": ["%s"], "sky_annulus": [450, 650], "object_diam": 300, "name": "%s", "plan": "A"},' % (file + "_SI.mat", name)
            print file, name, std

    out = out.rstrip(',')
    out += "]"

    f = open("stds.json", "w")

    f.write(out)
    f.close()
