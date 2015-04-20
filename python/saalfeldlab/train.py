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
    parser.add_argument( '--training-data', '-d', required=True, type=str, help='Path to training data (hdf5).' )
    parser.add_argument( '--data-type', '-t', default='f4', type=str, help='Datatype used for classifier.' )
    parser.add_argument( '--output', '-o', required=True, type=str, help='Path to file for storing classifier (hdf5).' )
    parser.add_argument( '--internal-path', '-i', required=True, type=str, help='Internal path for storing classifier.' )
    parser.add_argument( '--save-normalization', '-s', default='', type=str, help='If specified, data will be normalized and mean and variance will be save to file.' )
    parser.add_argument( '--save-normalization-internal', '-n', default='handler', type=str, help='Internal path for handler.' )
    parser.add_argument( '--subsample', '-b', default='-1', type=str, help='Number of data sample to use.' )
    parser.add_argument( '--kwargs', '-k', default='{}', type=str, help='Specify classifier parameters as dictionary.' )

    args = parser.parse_args()

    classifier = getattr( classifiers, args.classifier )( **ast.literal_eval( args.kwargs ) )

    handler      = data.DataHandler( args.training_data, np.dtype( args.data_type ), int( args.subsample ) )
    print 'reading the data'
    data, labels = handler.createFeatureMatrixAndLabelVector()
    print 'done'

    #if args.subsample >= 0 :
    #   print 'subsampling'
    #   data, labels = handler.subSampleData( data, labels )
    #   print 'done'

    if not args.save_normalization == '':
        data = handler.normalizeData( data )
        handler.save( args.save_normalization, args.save_normalization_internal )

    print 'training'
    classifier.train( data, labels )
    classifier.save( args.output, args.internal_path )
    print 'done'

    sys.exit( 0 )
