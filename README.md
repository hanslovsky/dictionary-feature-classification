dictionary-feature-classification
=================================

# Formatting data
This project accepts data in H5 format

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
