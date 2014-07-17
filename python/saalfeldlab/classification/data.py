import h5py
import numpy as np
import fnmatch
import os
import os.path


class DataHandler( object ):
    LABEL_PATH   = 'labels'
    SAMPLE_AXIS  = 0
    FEATURE_AXIS = 1
    
    def __init__(self, filename, dtype ):
        if ( os.path.isdir( filename )):
           self.dirname  = filename
           self.filename = "" 
        else:
           self.filename = filename

        self.dtype    = dtype
        self.mean     = None
        self.var      = None

    def createFeatureMatrixAndLabelVector( self ):
        if ( self.filename ):
			# we have a file, not a directory
            return self.createFeatureMatrixAndLabelVectorFile( 1 )
        else:
			# we have a directory of h5's
			self.filename = ""
			data   = np.array([]) 
			labels = np.array([])
			print "data is", repr(data)
			for f in os.listdir( self.dirname ):
				if fnmatch.fnmatch(f, '*.h5'):  
					self.filename = self.dirname + os.sep + f 
					thesedata,theselabels = self.createFeatureMatrixAndLabelVectorFile( 0 )
					if (data.size):	
						data = np.append( data, thesedata, axis=self.SAMPLE_AXIS)
						labels = np.append( labels, theselabels, axis=self.SAMPLE_AXIS)
					else:
						data   = thesedata;
						labels = theselabels;

			self.filename = ""
			return data, labels

    def createFeatureMatrixAndLabelVectorFile( self, enforce_labels = 1 ):
        with h5py.File( self.filename, 'r' ) as f:
            if not DataHandler.LABEL_PATH in f.keys():
                raise self.FileFormatException( "Label path %s not present in %s." % ( self.LABEL_PATH, self.filename ) )

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
                        raise self.DataInconsistencyException

                if ( enforce_labels and not (int(label) == idx + 1) ):
                    raise self.LabelInconsistencyException( "Labels not following 1,2,..." )

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
                #labelVector[ slicing[ DataHandler.SAMPLE_AXIS ], 0 ] = np.zeros( step ) + idx # int( label )
                labelVector[ slicing[ DataHandler.SAMPLE_AXIS ], 0 ] = np.zeros( step ) + int( label ) - 1 

                currentStart += step

            return featureMatrix, labelVector


    def normalizeData( self, data ):
        if self.mean is None:
            self.mean = np.mean( data, axis = DataHandler.SAMPLE_AXIS )
        if self.var is None:
            self.var  = np.var( data,  axis = DataHandler.SAMPLE_AXIS )

        # this works only if SAMPLE_AXIS == 0

        return ( data - self.mean ) / self.var

    def save( self, fileName, pathInFile ):
        with h5py.File( fileName, 'a' ) as f:
            if pathInFile in f.keys():
                del f[ pathInFile ]
            dg = f.create_group( pathInFile )
            for member in dir( self ):
                value = getattr( self, member )
                if callable( value ) or member.startswith( '__' ) or member[0].isupper() or member == 'dtype':
                    continue
                if value is None:
                    dg.create_dataset( name = member, shape=(1,) )
                    dg.attrs[ 'None' ] = [ 'True' ]
                else:
                    dg.create_dataset( name = member, data = value )



    def load( self, fileName, pathInFile ):
        with h5py.File( fileName, 'r' ) as f:
            dg = f[ pathInFile ]
            for member in dg.keys():
                if 'None' in dg[member].attrs.keys():
                    setattr( self, member, None )
                else:
                    setattr( self, member, dg[member].value )
                
                  


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
