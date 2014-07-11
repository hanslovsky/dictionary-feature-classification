import os
import sys
sys.path.append( os.path.dirname( os.path.realpath( __file__ ) ) )

from classification import *
import numpy as np
import vigra
import itertools
import inspect
import argparse
import ast


if __name__ == "__main__":
    availableClassifiers = [ c for c in dir( classifiers ) if inspect.isclass( getattr( classifiers, c ) ) and issubclass( getattr( classifiers, c ), classifiers.Classifier ) and not c == 'Classifier' ]
    parser = argparse.ArgumentParser()
    parser.add_argument( '--classifier', '-c', required=True, type=str, help='Specify classifier (choose from %s).' % availableClassifiers )
    parser.add_argument( '--classifier-file', '-f', required=True, type=str, help='Specify file from which classifier is loaded.' )
    parser.add_argument( '--classifier-path', '-p', required=True, type=str, help='Internal path for classifier.' )
    parser.add_argument( '--data', '-d', required=True, type=str, help='Path to training data (hdf5).' )
    parser.add_argument( '--data-type', '-t', default='f4', type=str, help='Datatype used for classifier.' )
    parser.add_argument( '--output', '-o', required=True, type=str, help='Path to file for storing prediction (hdf5).' )
    parser.add_argument( '--output-internal', '-i', default='prediction', help='Internal path for storing prediction.' )
    parser.add_argument( '--load-normalization', '-l', default='', type=str, help='If specified, data will be normalized with mean and variance loaded from file.' )
    parser.add_argument( '--load-normalization-internal', '-n', default='handler', type=str, help='Internal path for handler.' )
    parser.add_argument( '--predict-probabilities', '-P', action='store_true', help='Predict probabilities instead of labels.' )

    args = parser.parse_args()

    classifier = getattr( classifiers, args.classifier )( )
    classifier.load ( args.classifier_file, args.classifier_path )

    handler      = data.DataHandler( args.data, np.dtype( args.data_type ) )
    data, labels = handler.createFeatureMatrixAndLabelVector()

    if not args.load_normalization == '':
        handler.load( args.load_normalization, args.load_normalization_internal )
        data = handler.normalizeData( data )
        

    if args.predict_probabilities:
        prediction = classifier.predictProbabilities( data )
    else:
        prediction = classifier.predictLabels( data )

    correct = np.sum( prediction.flat == labels.flat )

    print correct, np.product( prediction.shape ), correct * 1.0 / np.product( prediction.shape )
        

    sys.exit( 0 )
