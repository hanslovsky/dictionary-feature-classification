import h5py
import numpy as np

class DataHandler( object ):
    LABEL_PATH   = 'labels'
    SAMPLE_AXIS  = 0
    FEATURE_AXIS = 1
    
    def __init__(self, filename, dtype ):
        self.filename = filename
        self.dtype    = dtype

    def createFeatureMatrixAndLabelVector( self ):
        with h5py.File( self.filename, 'r' ) as f:
            if not DataHandler.LABEL_PATH in f.keys():
                raise FileFormatException( "Label path %s not present in %s." % ( LABEL_PATH, self.filename ) )

            labelDir  = f[ DataHandler.LABEL_PATH ]
            keys      = sorted( labelDir.keys(), key = lambda x: int( x ) )
            nSamples  = 0
            nFeatures = None
            

            # check for consistency and get size of features
            for idx, label in enumerate( keys ):
                data = labelDir[ label ]
                if ( idx == 0 ):
                    nFeatures = data.shape[ DataHandler.FEATURE_AXIS ]
                else:
                    if not nFeatures == data.shape[ DataHandler.FEATURE_AXIS ]:
                        raise DataInconsistencyException

                if not int(label) == idx + 1:
                    raise LabelInconsistencyException( "Labels not following 1,2,..." )

                nSamples += data.shape[ DataHandler.SAMPLE_AXIS ]

                
            featureMatrixShape = [ 0, 0 ]
            
            featureMatrixShape[ DataHandler.SAMPLE_AXIS ]  = nSamples
            featureMatrixShape[ DataHandler.FEATURE_AXIS ] = nFeatures
            
            featureMatrix = np.empty( tuple( featureMatrixShape ), dtype=self.dtype )
            labelVector   = np.empty( ( nSamples, 1 ), dtype=np.uint32 )

            currentStart = 0

            slicing = [ slice( None ), slice( None ) ]
            
            for idx, label in enumerate( keys ):
                data = labelDir[label]
                step = data.shape[ DataHandler.SAMPLE_AXIS ]

                slicing[ DataHandler.SAMPLE_AXIS ] = slice( currentStart, currentStart + step )

                featureMatrix[ slicing ] = data[...]
                labelVector[ slicing[ DataHandler.SAMPLE_AXIS ], 0 ] = np.zeros( step ) + idx # int( label )

                currentStart += step

            return featureMatrix, labelVector
                
                  


    class DataInconsistencyException( Exception ):
        pass

    class LabelInconsistencyException( Exception ):
        pass

    class FileFormatException( Exception ):
        pass



if __name__ == "__main__":

    handler = DataHandler( "/nobackup/saalfeld/john/forPhilipp/testFeatures.h5", np.float32 )
    featureMatrix, labelVector = handler.createFeatureMatrixAndLabelVector()

    print featureMatrix.shape,'\n', labelVector.shape
