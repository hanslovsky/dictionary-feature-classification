import numpy as np
import vigra



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
