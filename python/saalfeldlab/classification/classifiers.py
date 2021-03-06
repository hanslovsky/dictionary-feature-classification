import numpy as np

import vigra

from sklearn import ensemble
from sklearn import naive_bayes
from sklearn import neighbors
from sklearn import svm

from sklearn.externals import joblib

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

    def save( self, filePath, pathInFile ):
        pass

    def load( self, filePath, pathInFile ):
        pass


class RandomForest( Classifier ):
    
    def __init__( self, **kwargs ):
        if 'randomSeed' in kwargs.keys():
            self.seed = kwargs[ 'randomSeed' ]
            del kwargs['randomSeed']
        else:
            self.seed = 10
        self.rf = vigra.learning.RandomForest( **kwargs )
        self.oob = None

    def train( self, data, labels ):
        self.oob = self.rf.learnRF( data, labels, int ( self.seed ) )
        return self

    def predictLabels( self, data ):
        return self.rf.predictLabels( data )

    def predictProbabilities( self, data ):
        return self.rf.predictProbabilities( data )

    def save( self, filePath, pathInFile ):
        self.rf.writeHDF5( filePath, pathInFile )

    def load( self, filePath, pathInFile ):
        self.rf = vigra.learning.RandomForest( filePath, pathInFile )


class SVM( Classifier ):

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

    def save( self, filePath, pathInFile ):
        joblib.dump( self.classifier, filePath )
    
    def load( self, filePath, pathInFile ):
        self.classifier = joblib.load( filePath )


class GaussianNaiveBayes( Classifier ):

    def __init__( self, **kwargs ):
        self.kwargs = kwargs
        self.classifier = naive_bayes.GaussianNB( **self.kwargs )

    def train( self, data, labels ):
        self.classifier.fit( data, labels.flat )
        return self

    def predictLabels( self, data ):
        return self.classifier.predict( data )[ :, np.newaxis ]

    def predictProbabilities( self, data ):
        return self.classifier.predict_proba( data )

    def save( self, filePath, pathInFile ):
        joblib.dump( self.classifier, filePath )

    def load( self, filePath, pathInFile ):
        self.classifier = joblib.load( filePath )


class KNearestNeighbor( Classifier ):

    def __init__( self, **kwargs ):
        self.kwargs = kwargs
        self.classifier = neighbors.KNeighborsClassifier( **self.kwargs )

    def train( self, data, labels ):
        self.classifier.fit( data, labels.flat )
        return self

    def predictLabels( self, data ):
        return self.classifier.predict( data )[ :, np.newaxis ]

    def predictProbabilities( self, data ):
        return self.classifier.predict_proba( data )

    def save( self, filePath, pathInFile ):
        joblib.dump( self.classifier, filePath )

    def load( self, filePath, pathInFile ):
        self.classifier = joblib.load( filePath )


class AdaBoost( Classifier ):

    def __init__( self, **kwargs ):
        self.kwargs = kwargs
        self.classifier = ensemble.AdaBoostClassifier( **self.kwargs )

    def train( self, data, labels ):
        self.classifier.fit( data, labels.flat )
        return self

    def predictLabels( self, data ):
        return self.classifier.predict( data )[ :, np.newaxis ]

    def predictProbabilities( self, data ):
        return self.classifier.predict_proba( data )

    def save( self, filePath, pathInFile ):
        joblib.dump( self.classifier, filePath )

    def load( self, filePath, pathInFile ):
        self.classifier = joblib.load( filePath )
