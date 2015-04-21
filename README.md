dictionary-feature-classification
=================================

# Classification 
The following trains a classifier on a training set, predicts classes for the training and test data, and
reports accuracy metrics.
```bash
$ classifyAndTest training_data.h5 test_data.h5 output_dir classifier_type subsample_param
```
See below for a details on the format of the expected h5 files.

## Subsample
One can subsample the input data during training to speed-up the process.
The default behavior uses all the data for training ( subSample = "-1" ).
Specifying a positive integer (N) will use that number of training samples.
Specifying a float value in (0,1] will use that fraction of training samples.

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
```bash
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
