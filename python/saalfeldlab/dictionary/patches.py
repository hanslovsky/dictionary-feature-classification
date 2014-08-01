import numpy as np
import itertools
import math

def getPatches( img, patchSize, N=-1, order='F'):

    pNumel = np.prod(patchSize)
    pRad   = (patchSize-1)/2
    nDim   = patchSize.shape[0]

    sz = img.shape 

    if N > 0 :

        startIdx = pRad 
        endIdx   = sz - pRad - 1
        coords = np.zeros( [nDim, N] )

        for d in range( 0, nDim ):
           coords[d,:] = np.random.random_integers( startIdx[d], endIdx[d], N ) 

        patchList = np.zeros( [N, pNumel], order=order)
        for i in range(0,N):
            thisPatch = img[ coords[0,i]-pRad[0] : coords[0,i]+pRad[0]+1, 
                             coords[1,i]-pRad[1] : coords[1,i]+pRad[1]+1, 
							 coords[2,i]-pRad[2] : coords[2,i]+pRad[2]+1 ]

            patchList[i,:] = thisPatch.flatten()
    else:
        # find the number of elements
        subSz = np.array(sz) - 2 * pRad
        N = np.prod( subSz ) 
        patchList = np.zeros( [N, pNumel], order=order)
        print " subSz + pRad ", subSz + pRad
        for i,(x,y,z) in enumerate( itertools.product( *map( xrange, subSz ))):
            print i,x,y,z
            x+=pRad[0]
            y+=pRad[1]
            z+=pRad[2]
            print x
            print y
            print z
            print x+patchSize[0]
            print y+patchSize[1]
            print z+patchSize[2]
            thisPatch = img[ x : x+patchSize[0], 
                             y : y+patchSize[1], 
							 z : z+patchSize[2] ]

            print thisPatch.shape

            patchList[i,:] = thisPatch.flatten()

    return patchList

def jointEntropy( hiPatch, loPatch, downsampleFactor, bitDepth ):
    """ Joint entropy 
    """


def conditionalEntropy( loPatch, downsampleFactor, bitDepth ):
    """ Entropy (H)
    """ 
    lut,pmf = generateLowToHiResLutPmf( downsampleFactor, bitDepth )
   
    H = 0.0 # the conditional entropy

    for value in loPatch.flatten():
        print value
        H += pmf[ value ] * math.log( pmf[ value ] )


def patchHypotheses( loPatch, downsampleFactor ):
    """ Expect loPatch to be 3d 
    """
    print "gen patch hypotheses"
    
    
def generateLowToHiResLutPmf( D, bitDepth ):
    """ Generate a look-up-table that gives all possible high res
    configurations that could give rise to a given low-res value.

        For every level at the input bitDepth v, generates a list of 
        D values also at that bitDepth whose average is v.
    """
    print "generateLowToHiResLUT for ", D, " downsample factor"

    numLevels = 2**bitDepth
    print "num levels: ", numLevels    

    # Generate all D*numLevels combinations of values
    s = np.indices( tuple( numLevels * np.ones(D,)))
    t = np.reshape( s, ( s.shape[0], int(np.prod( s.shape[1:]))))
   
    # Get the lo res values
    lut = np.round( np.mean( t, axis=0 ))
    pmf = np.float32(np.bincount( np.int64(lut) )) / lut.shape[0]

    return lut,pmf

    
