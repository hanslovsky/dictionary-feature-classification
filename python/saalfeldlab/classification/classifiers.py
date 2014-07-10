import numpy as np

import vigra

from sklearn import svm

import h5py



class Classifier( object ):

    def __init__( self ):
        raise Exception( "Not implemented!" )

    def train( self, data, labels ):
        pass

    def predictLabels( self, data ):
        pass

    def predictProbabilities( self, data ):
        pass

    def writeToHdf5( self, filePath, pathInFile ):
        pass

    def readFromHdf5( self, filePath, pathInFile ):
        pass


class RandomForest( Classifier ):
    
    def __init__( self, randomSeed, **kwargs ):
        self.rf = vigra.learning.RandomForest( **kwargs )
        self.seed = randomSeed
        self.oob = None

    def train( self, data, labels ):
        print data.dtype, data.shape, labels.dtype, labels.shape
        self.oob = self.rf.learnRF( data, labels ) #, int ( self.seed ) )
        return self

    def predictLabels( self, data ):
        return self.rf.predictLabels( data )

    def predictProbabilities( self, data ):
        return self.rf.predictProbabilities( data )

    def writeToHdf5( self, filePath, pathInFile ):
        self.rf.writeHDF5( filePath, pathInFile )

    def readFromHdf5( self, filePath, pathInFile ):
        self.rf = vigra.learning.RandomForest( filePath, pathInFile )



class SVM( object ):

    def __init__( self, **kwargs ):
        self.kwargs = kwargs
        self.classifier = svm.SVC()

    def train( self, data, labels ):
        self.classifier.fit( data, labels.flat )
        return self

    def predictLabels( self, data ):
        return self.classifier.predict( data )[ :, np.newaxis ]

    def predictProbabilities( self, data ):
        return self.predictLabels( data )

    def writeToHdf5( self, filePath, pathInFile ):
        with h5py.File( filePath, 'w-' ) as f:
            group  = f.create_group( pathInFile )
            params = self.classifier.get_params( deep = False )
            for p, v in params.iteritems():
                group.create_dataset( p, data = np.array( v ) )

    def readFromHdf5( self, filePath, pathInFile ):
        params = {}
        with h5py.File( FilePath, 'r' ) as f:
            group = f[ pathInFile ]
            for p in group.keys():
                params[p] = group[p][...]

        self.classifier.set_params( **params )
