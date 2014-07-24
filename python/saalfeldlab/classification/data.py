import h5py
import numpy as np
import fnmatch
import os
import os.path


class DataHandler( object ):
    LABEL_PATH   = 'labels'
    SAMPLE_AXIS  = 0
    FEATURE_AXIS = 1
    
    def __init__(self, filename, dtype, sampleParam=-1 ):
        if ( os.path.isdir( filename )):
           self.dirname  = filename
           self.fList    = os.listdir( self.dirname )
           self.filename = "" 
        else:
           self.filename = filename

        self.dtype    = dtype
        self.mean     = None
        self.var      = None

        self.doSubsample = False

        # params and helper objects for random sampling
        if (sampleParam > 0 and sampleParam <=1 ):
            self.sampleRatio = sampleParam
            self.doSubsample = True
        elif (sampleParam > 1):
            self.N           = int(sampleParam)
            self.sampleRatio = None
            self.doSubsample = True
        else:
            self.sampleRatio = None

        if ( self.doSubsample ):
		   self.randomSamplingSetup()


    def sampleCountFile( self, filename ):
        """ Count the number of samples per class in an hdf5 data file. 

        Returns a numpy array where the ith entry gives the number of 
        samples with label i in the file.
        """

        sampleCounts = np.zeros( 2, 'int64' ) 
        # we have a file, not a directory
        with h5py.File( filename, 'r' ) as f:
            #print f
            if not DataHandler.LABEL_PATH in f.keys():
                raise self.FileFormatException( "Label path %s not present in %s." % ( self.LABEL_PATH, filename ) )
            labelDir  = f[ DataHandler.LABEL_PATH ]
            keys      = sorted( labelDir.keys(), key = lambda x: int( x ) )
            
            # check for consistency and get size of features
            for idx, label in enumerate( keys ):
                shp = f[ self.LABEL_PATH + '/' + str(label)].shape
                sampleCounts[idx] += int( shp[self.SAMPLE_AXIS] )

        return sampleCounts

    def randomSamplingSetup( self ):
        """ Do prerequisite computations for random sampling of data. 
        """

        if ( self.filename ):
			# we have a file, not a directory
            self.sampleCounts = self.sampleCountFile( self.filename )
        else:
            self.sampleCounts = np.zeros([len(self.fList), 2], 'int64')
            i = 0
            #TODO make work for arbitrary label count    
            for f in self.fList:
                if fnmatch.fnmatch(f, '*.h5'):  
                    filename = self.dirname + os.sep + f 
                    self.sampleCounts[i,:] = self.sampleCountFile( filename )
                    i=i+1

            # cumulative sample counts 
            self.cumSampleCounts = np.cumsum(self.sampleCounts,0)
       
        if ( self.sampleRatio ):
           self.N = np.round( self.sampleRatio * self.cumSampleCounts[-1,:] )
        elif not self.N is np.ndarray:
           self.N = np.int_(np.round(self.N * (np.float_(self.cumSampleCounts[-1,:])/np.sum(self.cumSampleCounts[-1,:]))))


        i = np.random.permutation( int(self.cumSampleCounts[-1,0] ))
        j = np.random.permutation( int(self.cumSampleCounts[-1,1] ))
        i = np.sort( i[0:int(self.N[0])] )
        j = np.sort( j[0:int(self.N[1])] )

        self.randomSamples = [i,j]
        
        n = self.sampleToFileIdx( i, 0 )
        m = self.sampleToFileIdx( j, 1 )

        self.fileToIdxs = [n,m]

    def sampleToFileIdx( self, i, label):
        """ Requires that randomSamplingSetup be run """
        res = np.zeros( i.shape )
        for n in range(0, len(i)):
            res[n] = np.where( self.cumSampleCounts[:,label] > i[n] )[0][0]

        return res

    def fileIdxSamples( self, thisFileIndex, label ):
        j = np.where( self.fileToIdxs[label] == thisFileIndex )[0]
        res = self.randomSamples[label][j] 
        if thisFileIndex > 0 :
            res -= self.cumSampleCounts[ thisFileIndex-1, label ] 
        
        return res 

    def subSampleData( self, data, labels ):
        totSamples = data.shape[ self.SAMPLE_AXIS ]
        i = np.random.permutation( totSamples )
        i = i[ 0:self.N ]
        data = data[ i, : ]
        labels = labels[ i ]

        return data, labels

    def createFeatureMatrixAndLabelVector( self ):
        if ( self.filename ):
			# we have a file, not a directory
            return self.createFeatureMatrixAndLabelVectorFile( 1 )
        else:
            # we have a directory of h5's	
            self.filename = ""	
            data   = np.array([]) 
            labels = np.array([])
            #print "data is", repr(data)
            for i,f in enumerate( self.fList ):
                #print f
                if fnmatch.fnmatch(f, '*.h5'):  
                    self.filename = self.dirname + os.sep + f 
                    if self.doSubsample:
                       thesedata,theselabels = self.createFeatureMatrixAndLabelVectorFile( 0, i )
                    else:
                       thesedata,theselabels = self.createFeatureMatrixAndLabelVectorFile( 0 )

                    if thesedata is None: 
                        continue

                    if (np.prod(data.shape)):	
                        data = np.append( data, thesedata, axis=self.SAMPLE_AXIS)
                        labels = np.append( labels, theselabels, axis=self.SAMPLE_AXIS)
                    else:
				        data   = thesedata;
				        labels = theselabels;

                    #print data.shape

            self.filename = ""
            return data, labels


    # TODO implement random sampling for single file
    def createFeatureMatrixAndLabelVectorFile( self, enforce_labels = 1, file_index = -1 ):
        with h5py.File( self.filename, 'r' ) as f:
            if not DataHandler.LABEL_PATH in f.keys():
                raise self.FileFormatException( "Label path %s not present in %s." % ( self.LABEL_PATH, self.filename ) )

            labelDir  = f[ DataHandler.LABEL_PATH ]
            keys      = sorted( labelDir.keys(), key = lambda x: int( x ) )
            nSamples  = 0
            nFeatures = None
            
            # check for consistency and get size of features
            dataRowList = []
            for idx, label in enumerate( keys ):
                data = labelDir[ label ]
                if ( idx == 0 ):
                    nFeatures = data.shape[ DataHandler.FEATURE_AXIS ]
                else:
                    if not nFeatures == data.shape[ DataHandler.FEATURE_AXIS ]:
                        raise self.DataInconsistencyException

                if ( enforce_labels and not (int(label) == idx + 1) ):
                    raise self.LabelInconsistencyException( "Labels not following 1,2,..." )

                if self.doSubsample:
                   dataRows = self.fileIdxSamples( file_index, idx )
                   nSamples += len( dataRows )
                   if len(dataRows)==0:
                      dataRowList.append( np.array([-1]) )
                   else:
                      dataRowList.append( dataRows )

                   #print 'dataRowList: ', dataRowList
                else:
                   nSamples += data.shape[ DataHandler.SAMPLE_AXIS ]

            if not nSamples:
               return None, None

            featureMatrixShape = [ 0, 0 ]
            featureMatrixShape[ DataHandler.SAMPLE_AXIS ]  = nSamples
            featureMatrixShape[ DataHandler.FEATURE_AXIS ] = nFeatures
            
            featureMatrix = np.empty( tuple( featureMatrixShape ), dtype=self.dtype )
            labelVector   = np.empty( ( nSamples, 1 ), dtype=np.uint32 )

            currentStart = 0
            slicing = [ slice( None ), slice( None ) ]

            #print 'nSamples: ', nSamples
            #print 'dataRowList: ', dataRowList
            
            for idx, label in enumerate( keys ):
                if self.doSubsample and dataRowList[idx][0] == -1 :
                    continue

                data = labelDir[label]

                if self.doSubsample:
                    step = len(dataRowList[idx])
                else:
                    step = data.shape[ DataHandler.SAMPLE_AXIS ]

                slicing[ DataHandler.SAMPLE_AXIS ] = slice( currentStart, currentStart + step )

                if self.doSubsample:
                    if file_index == -1:
                        raise self.ParameterException

                    # TODO return if dataRows is empty
                    featureMatrix[ slicing ] = data[ dataRowList[idx], : ]
                else:
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

    class ParameterException( Exception ):
        pass


if __name__ == "__main__":

    handler = DataHandler( "/nobackup/saalfeld/john/forPhilipp/testFeatures.h5", np.float32 )
    featureMatrix, labelVector = handler.createFeatureMatrixAndLabelVector()

    print featureMatrix.shape,'\n', labelVector.shape
