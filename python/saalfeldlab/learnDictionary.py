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
    #parser.add_argument( '--test-image', '-s', default="", type=str, 
    #      help='The test hdf5 image' )
    #parser.add_argument( '--downsample-factor', '-f', default=3, type=int, 
    #      help='Downsampling factor in z' )
    parser.add_argument( '--output', '-o', default="", type=str, 
          help='The output hdf5 file where the dictionary will be stored' )
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
#    parser.add_argument( '--objective-function', '-j', default=False, type=bool,
#          help="Compute and print value of objective function")
    parser.add_argument( '--threads', '-r', default=1, type=int, 
          help='Number of threads ' )
    parser.add_argument( '--verbose', '-v', dest='verbose', action='store_true',
          help='Verbose output' )

    args = parser.parse_args()
    K = args.dictionary_size
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

    params = {          'K' : K,
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


    # Write the dictionary
    if args.output:
        print "Writing dictionary"

        h5out = h5py.File( args.output, 'a' )
        h5out.create_dataset("dict", data=D)
        h5out.create_dataset("objective", data=R)

        # save parameters
        paramGroup = h5out.create_group("param")
        for pKey in params.keys():
            paramGroup.create_dataset( pKey, data=params[pKey])

        # save image and number of samples
        paramGroup.create_dataset( 'source_img', data=args.image )
        paramGroup.create_dataset( 'numSamples', data=N )

        h5out.flush()
        h5out.close()
    
        
    sys.exit( 0 )
