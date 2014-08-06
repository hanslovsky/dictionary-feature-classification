import os
import sys
import time
sys.path.append( os.path.dirname( os.path.realpath( __file__ ) ) )

from dictionary import patches
import numpy as np
import h5py
import spams

import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--image', '-i', required=True, type=str, 
          help='The source hdf5 image' )
    parser.add_argument( '--dictionary', '-c', default="", type=str, 
          help='The dictionary/codebook hdf5.' )
    parser.add_argument( '--image-internal-directory', '-d', default='main', type=str, 
          help='Internal directory of hdf5 image' )
    parser.add_argument( '--number', '-n', default='-1', type=int, 
          help='Number of patches' )
    parser.add_argument( '--patch-size', '-p', default='5-5-5', type=str, 
          help='Patch size' )
    parser.add_argument( '--threads', '-r', default=1, type=int, 
          help='Number of threads ' )
    parser.add_argument( '--verbose', '-v', dest='verbose', action='store_true',
          help='Verbose output' )

    args = parser.parse_args()
    N = args.number

    patchSize = np.fromstring( args.patch_size, dtype=int, sep='-')
    print "patch size: ", patchSize

    print "reading data"
    print args.image , "/", args.image_internal_directory

    f = h5py.File( args.image )
    im = f[ args.image_internal_directory ][...]

    print "sampling ", N, " patches"
    tic = time.time()
    X = patches.getPatches( im, patchSize, N )
    toc = time.time()
    t = toc - tic
    print 'time to grab patches: %f' % t 

    print X.shape
    print np.isfortran( X )

    params = {    'lambda1' : 0.15, 
               'numThreads' : args.threads,
                  'verbose' : args.verbose }

                  
    print "loading dictionary"
    tic = time.time()
    dfn = h5py.File( args.dictionary )	
    D =  dfn[ 'dict' ][...]
    D = np.asfortranarray( D )
    toc = time.time()	
    t = toc - tic
    print 'time to load dictionary: %f' % t 

    print D.shape
    print np.isfortran( D )

    #if args.objective_function:
    if True:
        print "computing objective function"

        ###############################
        lst = [ 'L','lambda1','lambda2','mode','pos','ols','numThreads','length_path','verbose','cholesky']
        lparam = {'return_reg_path' : False}
        for x in lst:
            if x in params:
                lparam[x] = params[x]
        ###############################

        tic = time.time()

        alpha = spams.lasso( X, D = D, **lparam )
        xd = X - D * alpha 
        R = np.mean(0.5 * (xd * xd).sum(axis=0) + params['lambda1'] * np.abs(alpha).sum(axis=0))
        toc = time.time()

        print "  objective function value: %f" % R
        t = toc - tic
        print 'time of computation for objective function: %f' % t
  

    sys.exit( 0 )
