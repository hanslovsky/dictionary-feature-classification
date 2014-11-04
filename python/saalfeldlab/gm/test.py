import numpy as np

import opengm

import sets

import scipy.spatial
import scipy.io as io

import vigra


def create_gm( unaries, regularizer ):
    gm = opengm.grid2d2Order( unaries, regularizer )


def dictionaryToList( dictionary ):
    result = []
    for index in xrange( dictionary.shape[ 1 ] ):
        result.append( dictionary[ :, index ] )

    return result


def findNBest( patch3D, dict2D, N, shape, direction ):
    N       = min( N, len(dict2D ) )
    shape2D = ()
    
    for i in xrange(len(shape)):
        if i == direction:
            continue
        shape2D = shape2D + ( shape[i], )
        
    patch   = patch3D.reshape( shape )
    tree    = scipy.spatial.KDTree( dict2D, leafsize=1000 )
    slicing = [ slice(None) ] * len(shape)
    distances = np.empty( ( shape[direction], N ), dtype=np.float64 )
    indices   = np.empty( ( shape[direction], N ), dtype=np.uint32  )
    
    for index in xrange( shape[direction] ):
        slicing[ direction ] = index
        currentSlice = patch[ slicing ]
        currentSlice = currentSlice.reshape( 1, currentSlice.size )
        query = tree.query(currentSlice, k=N )
        distances[ index, ... ] = query[0]
        indices[ index, ... ]   = query[1]

    return distances, indices


def downsample( dict2D_List, shape, direction ):
    nPatches = len( dict2D_list )
    result   = np.empty( ( nPatches, shape[direction] * shape[2] ) )
    ratio    = shape[direction] * 1.0 / shape[2]
    slicing  = [ slice(None) ] * 2
    slicing[ direction ] = slice( None, None, ratio )
    
    sigma            = [ 0.0001, 0.0001 ]
    sigma[direction] = np.sqrt( ratio )
    
    for n in xrange( nPatches ):
        patch     = dict2D_list[ n ].reshape( shape[:2] )
        smoothed  = vigra.filters.gaussianSmoothing( patch, sigma=sigma )
        result[n,...] = smoothed[ slicing ].flatten()

    return result


def calculateDistance( i1, i2, dictionary, patches ):
    s = tuple( sets.Set( ( i1, i2 ) ) )
    if dictionary.has_key( s ):
        return dictionary[s]
    else:
        distance = np.sum( ( patches[i1] - patches[i2])**2 )
        dictionary[s] = distance
        return distance


def calculateDistances( indices1, indices2, dictionary, patches, alpha=1.0 ):
    result = np.zeros( ( indices1.shape[0], indices2.shape[0] ) )
    for e1, i1 in enumerate( indices1 ):
        for e2, i2 in enumerate( indices2 ):
            result[e1, e2] = alpha*calculateDistance( i1, i2, dictionary, patches )
    return result


def generateHighResPatch( state, indices, patches, direction, shape ):
    shapeHighRes  = tuple( shape[:2] ) + (state.shape[0],)
    result = np.empty( shapeHighRes )
    slicing = [slice(None)] * len(shape)
    for i, s in enumerate(state):
        slicing[direction] = i
        result[slicing] = patches[indices[i][s]].reshape(shape[:2])
    return result


if __name__ == "__main__":
    #dictionaries: /groups/saalfeld/home/bogovicj/projects/dictionary/forPhillip/
    #default values: http://stackoverflow.com/questions/12151306/argparse-way-to-include-default-values-in-help
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( "--two", '-t', help="Path to 2D dictionary.", required=True )
    parser.add_argument( "--three", '-T', help ="Path to 3D dictionary.", required=True )
    parser.add_argument( "--shape", '-s', help="Shape of 3D patches in the format x,y,z", required=True )
    parser.add_argument( "--direction", '-d', help="Axis along which patches are being matched.", required=True, type=int )
    parser.add_argument( "--number-of-neighbors", '-N', help='Number of neighbors for nearest neighbor search when generating unaries.', required=True, type=int )
    parser.add_argument( '--alpha', '-a', default=1.0, help='Weight of binary terms vs neighbor terms (default: %(default)s).', type=float )
    parser.add_argument( '--high-res', '-H', default='hr', help='Output directory for high-res dictionaries (default: %(default)s).', type=str )
    parser.add_argument( '--low-res',  '-l', default='lr', help='Output directory for low-res dictionaries (default: %(default)s).', type=str )

    args = parser.parse_args( )

    shape = [ int(x) for x in args.shape.split(',') ]

    dict2D = np.require( io.loadmat( args.two )['D2d'], dtype = np.float32 )
    dict3D = np.require( io.loadmat( args.three )['D'], dtype = np.float32 )

    dict2D_list = dictionaryToList( dict2D )
    dict3D_list = dictionaryToList( dict3D )

    dict2D_downsampled = downsample( dict2D_list, shape, args.direction )
    nPatches3D = len( dict3D_list )
    distances2D = {}

    nNeighbors = args.number_of_neighbors
    numLabels  = [ nNeighbors ] * shape[args.direction]

    highResPatches = []

    for p in xrange( nPatches3D ):

        patch = dict3D_list[p]
        distances, indices = findNBest( patch, dict2D_downsampled, nNeighbors, shape, args.direction )

        distances3D = distances[ np.newaxis, ... ]

        gm = opengm.graphicalModel( numLabels, operator='adder' )# opengm.grid2d2Order( distances3D, regularizer=Regularizer() )

        for v in xrange( distances.shape[ 0 ] ):
            functionId = gm.addFunction( distances[ v ] )
            gm.addFactor( functionId, v )
            if v > 0:
                binary = calculateDistances( indices[v-1], indices[v], distances2D, dict2D_list, args.alpha )
                binaryFunctionId = gm.addFunction( binary )
                gm.addFactor( binaryFunctionId, [v-1, v] )
                
        print 'before creating inference'
        # inference = opengm.inference.TreeReweightedBp( gm )
        inference = opengm.inference.DynamicProgramming( gm )
        print 'before inferring'
        inference.infer()
        print 'done inferring at time', p
        result = inference.arg()
        patchHR  = generateHighResPatch( result, indices, dict2D_list, args.direction, shape )
        vigra.impex.writeHDF5( patch.reshape( shape ), 'lr/%03d.h5' % p, 'data' )
        vigra.impex.writeHDF5( patchHR, 'hr/%03d.h5' % p, 'data' )
        

    

