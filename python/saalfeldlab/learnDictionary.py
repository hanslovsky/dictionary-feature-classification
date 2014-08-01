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
    parser.add_argument( '--image-internal-directory', '-d', default='main', type=str, 
          help='Internal directory of hdf5 image' )
    parser.add_argument( '--number', '-n', default='-1', type=int, 
          help='Number of patches' )
    parser.add_argument( '--dictionary-size', '-k', default=100,
          type=int, help='Dictionary size' )
    parser.add_argument( '--patch-size', '-p', default='5 5 5', type=str, 
          help='Patch size' )
    parser.add_argument( '--iterations', '-t', default=1000, type=int,
          help='Dictionary learning iterations' )
    parser.add_argument( '--batch-size', '-b', default=100, type=int,
          help="Batch size")
    parser.add_argument( '--threads', '-r', default=1, type=int, 
          help='Number of threads ' )
    parser.add_argument( '--verbose', '-v', default=False, type=bool, 
          help='Verbose output' )

    args = parser.parse_args()

    patchSize = np.fromstring( args.patch_size, dtype=int, sep=' ')
    print patchSize

    print "reading data"
    print args.image , "/", args.image_internal_directory

    f = h5py.File( args.image )
    im = f[ args.image_internal_directory ][...]

    M = args.dictionary_size
    N = args.number
    print "sampling ", N, " patches"
    X = patches.getPatches( im, patchSize, N )
    print X.shape

    params = {          'K' : M,
                  'lambda1' : 0.15, 
               'numThreads' : args.threads,
                'batchsize' : args.batch_size,
                     'iter' : args.iterations,
                  'verbose' : args.verbose }


    print "learning dictionary"
    tic = time.time()
    D = spams.trainDL( X, **params )
    
    toc = time.time()

    t = toc - tic
    print 'time of computation for Dictionary Learning: %f' % t 


