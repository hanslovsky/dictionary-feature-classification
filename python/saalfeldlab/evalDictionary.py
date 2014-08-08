import os
import sys
import time
sys.path.append( os.path.dirname( os.path.realpath( __file__ ) ) )

from dictionary import patches
from dictionary import evaluation
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
    parser.add_argument( '--do-upsampling', '-u', dest='do_upsampling', action='store_true',
          help='Do upsampling test' )
    parser.add_argument( '--downsample-factor', '-f', default=3, type=int, 
          help='Downsampling factor in z' )
    parser.add_argument( '--verbose', '-v', dest='verbose', action='store_true',
          help='Verbose output' )

    args = parser.parse_args()
    N = args.number
    doUpsamp = args.do_upsampling

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

    params = {    'lambda1' : 0.15, 
               'numThreads' : args.threads,
                  'verbose' : args.verbose }

    ###############################
    lst = [ 'L','lambda1','lambda2','mode','pos','ols','numThreads','length_path','verbose','cholesky']
    lparam = {'return_reg_path' : False}
    for x in lst:
        if x in params:
            lparam[x] = params[x]
    ###############################
                  
    print "loading dictionary"
    tic = time.time()
    dfn = h5py.File( args.dictionary )	
    D =  dfn[ 'dict' ][...]
    D = np.asfortranarray( D )
    # Read dictionary parameters
    #dparamsGrp = dfn[ 'params ']
    #tmp = []
    #for k in dparamsGrp.keys():
    #    tmp.append( (k, dparamsGrp[k]))    

    #lparam = dict( tmp )
    toc = time.time()	
    t = toc - tic
    print 'time to load dictionary: %f' % t 

    print "D shape",D.shape
    print "X shape",X.shape

    #if args.objective_function:
    #if True:
    if False:
        print "computing objective function"
        
        tic = time.time()

        alpha = spams.lasso( X, D = D, **lparam )
        xd = X - D * alpha 
        R = np.mean(0.5 * (xd * xd).sum(axis=0) + params['lambda1'] * np.abs(alpha).sum(axis=0))
        toc = time.time()

        print "  objective function value: %f" % R
        t = toc - tic
        print 'time of computation for objective function: %f' % t
 
    if doUpsamp:
        print "doing upsampling evaluation"
        tic = time.time()
        # dsz - size of downsampled image patch
        X_ds,dsz = evaluation.downsamplePatchList( X, patchSize, np.array([1,1,args.downsample_factor]))
        X_ds = np.asfortranarray( X_ds )

        D_ds,Ddsz = evaluation.downsamplePatchList( D, patchSize, np.array([1,1,args.downsample_factor]))  
        D_ds = np.asfortranarray( D_ds )
    
        print "dz patch sz ", dsz, " Ddsz ", Ddsz
        print "D_ds shape",D_ds.shape
        print "X_ds shape",X_ds.shape
        alpha = spams.lasso( X_ds, D = D_ds, **lparam )

        xd_ds = X - D * alpha
        R_ds = np.mean((xd_ds * xd_ds).sum(axis=0))
        toc = time.time()

        print "mean squared error: %f" % R_ds
        print 'time of computation: %f' % t


        print "comparing up NN upsampling"
		# Compare to linearly upsampled patch
        tic = time.time()
        X_us = evaluation.upsamplePatchList( X_ds, dsz, args.downsample_factor )
        xd_us = X - X_us
        R_us = np.mean( (xd_us * xd_us).sum(axis=0))
        toc = time.time()

        print "mean squared error upsample nearest: %f" % R_us
        print 'time of computation: %f' % t


    sys.exit( 0 )
