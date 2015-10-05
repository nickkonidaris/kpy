
import argparse, os, pdb, sys
import numpy as np
import pyfits as pf



def roll_file(infile, X, Y):
    ''' Load fits infile and roll it in X,Y '''
    
    X = np.int(np.round(X))
    Y = np.int(np.round(Y))

    FF = pf.open(infile)
    FF[0].data = np.roll(FF[0].data, X, axis=1)
    FF[0].data = np.roll(FF[0].data, Y, axis=0)

    FF[0].header['ROLLX'] = (X, "Rolled in X")
    FF[0].header['ROLLY'] = (Y, "Rolled in Y")


    FF.writeto(infile, clobber=True)

    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Roller.py:

            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('infile', type=str, help='File to roll')
    parser.add_argument('--X', type=float, default=0, help='Roll amount X')
    parser.add_argument('--Y', type=float, default=0, help='Roll amount Y')
    args = parser.parse_args()

    
    roll_file(args.infile, args.X, args.Y)


