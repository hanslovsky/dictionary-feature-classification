import numpy as np
from sklearn.metrics import *

def evalProb( truth, prob, writeOutput = True, suffix="" ):
    precision, recall, thresholds = precision_recall_curve( truth, prob )
    fMeasure = ( 2 * precision * recall ) / ( precision + recall )
    auPRcurve = auc( recall, precision )
    print " ################################# "

    # Save the curves to files
    if( writeOutput ):
        np.savetxt( "precision_"+suffix+".csv", precision, delimiter=',', fmt='%0.4f' )
        np.savetxt( "recall_"+suffix+".csv", recall, delimiter=',', fmt='%0.4f' )
        np.savetxt( "thresholds_"+suffix+".csv", thresholds, delimiter=',', fmt='%0.4f' )
        np.savetxt( "fMeasure_"+suffix+".csv", fMeasure, delimiter=',', fmt='%0.4f' )

    # find the best fMeasure
    bestF_index = np.argmax( fMeasure )
    bestF = fMeasure[ bestF_index ]
    bestPrecision = precision[ bestF_index ]
    bestRecall = recall[ bestF_index ]
    bestThreshold = recall[ bestF_index ]
    
    print "     F\t Prec\t Recl\tThresh"
    print "%0.4f\t%0.4f\t%0.4f\t%0.4f" % ( bestF,  bestPrecision, bestRecall, bestThreshold )
    print " ################################# "
    print "AUC: %0.4f" % ( auPRcurve )

    # evaluate the binary measures using the bestThreshold
    evalClass( truth, ( prob > bestThreshold ))


def evalClass( truth, pred ):
    print( classification_report( truth, pred ))

    correct = np.sum( pred.flat == truth.flat )
    acc = correct * 1.0 / np.product( pred.shape )

    labelList = np.unique( truth ) 
    # Compute balanced class accuracy
    acc_bal = 0.0
    for l in labelList:
       numL = np.sum( truth.flat == l )
       # the true positive rate for this label
       acc_bal += np.sum( np.logical_and( (truth.flat == l), (pred.flat == l))) 

    # divide by number of labels
    acc_bal /= np.prod( truth.shape )

    print "Ncorrect\tN\tacc\tacc_bal"
    print "%d\t\t%d\t%0.4f\t%0.4f" % ( correct,  np.product( pred.shape ), acc, acc_bal ) 
