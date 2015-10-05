
import shutil
import glob
import os

def go(fr, to):
    '''Current as of Jan 25 2014, move raw files to proper directory structure'''
    mnths = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", 
        "oct", "nov", "dec"]
    
    if type(fr) == str:
        files = glob.glob(fr)
    else:
        files = fr

    for path in files:
        fname = path.split("/")[-1]
        try:
            td, h, min, s = fname.split("_")

            if td[0:3] == 'ifu':
                otype = 'ifu'
                y,m,d = td[3:7], td[7:9], td[9:11]
            else:
                otype = 'rc'
                continue
                y,m,d = td[2:6], td[6:8], td[8:10]

        except:
            print "Skipping %s" % fname
            continue

        
        if int(h) >= 12: outday = d
        else: outday = "%2.2i" % (int(d)-1)

        outdir = os.path.join(to, y+ mnths[int(m)-1]+outday)
        outfile = otype + y+m+d+"_"+h+"_"+min+"_"+s

        if os.path.exists(outdir) == False:
            os.mkdir(outdir)
        print path, outfile
        try: shutil.copy2(path, os.path.join(outdir,outfile))
        except IOError: continue

        
if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        raise Exception("not enough arguments")
        

    go(sys.argv[1:-1], sys.argv[-1])
        

