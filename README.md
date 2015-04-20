dictionary-feature-classification
=================================

# Formatting data
This project accepts data in the [HDF5 format](https://hdfgroup.org/HDF5/)

## Single file 
* The H5 file must have a 'labels' group with one Dataset per class.
* Currently classes must be named  "1","2",..."N"
* The Dataset for each class must be 2d ( N x M )
  * N - number of observations
  * M - number of variables

### Example
The below is an example of the format for the famous [iris data](http://archive.ics.uci.edu/ml/datasets/Iris).
```
$ h5ls iris_all.h5
labels                   Group
$ h5ls iris_all.h5/labels
1                        Dataset {50, 4}
2                        Dataset {50, 4}
3                        Dataset {50, 4}
```

## Multiple files
It is sometimes convenient to store data in multiple h5 files.  In this case, 
one can specify a directory containing multiple h5 files.  In this case, all h5
files in that directory will be included in the dataset.  
