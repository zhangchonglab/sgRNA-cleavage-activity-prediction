# this script is used to train a gradient boosting regression model for sgRNA activity prediction

from __future__ import print_function
import os
import sys
import numpy as np
import pandas as pd
import logging
import argparse
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingRegressor  #GBM algorithm
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pickle
import math
plt.rcParams["font.family"] = "Times New Roman"

def GBR_pred_parseargs():
    """
    Parse arguments. Only used when GBR_pred.py is executed directly.
    """
    parser=argparse.ArgumentParser(description='Gradient Boosting Regression to deciper sgRNA activity from sequence feature')
    parser.add_argument('-d', '--dataset', required=True, help='raw dataset of sgRNA ID and feature')
    parser.add_argument('-i', '--id', required=True, default='sgRNAID', help='column of sgRNA ID')
    parser.add_argument('-m', '--model', required=True, default='saved_model/Cas9_sgRNA_activity_GBR.pickle', help='pickle GBR model')
    parser.add_argument('-n', '--normalization', required=True, default='Yes', help='whether to perform normalization on results or not')
    parser.add_argument('-p', '--prefix', required=True, help='prefix of the result file')
    args=parser.parse_args()
    return args

# ///////////////////////////////////////////////////////
def perform_normalization(activity):
    """
    linear project values in the activity array to [0, 1]
    Parameters:
    ======================
    activity: array of scores
    """
    max_value=np.nanmax(activity)
    min_value=np.nanmin(activity)
    scaler=max_value-min_value
    return (activity-min_value)/scaler

# //////////////////////////////////////////
def GBR_pred_main(args):
    """
    Main entry
    """
    """
    get input parameters
    """
    # raw dataset with features and scores
    my_data=pd.read_csv(args.dataset)
    # id column in this dataset
    id_column=args.id
    # trained models used for sgRNA activity prediction
    model=args.model
    gbm_tuned=pickle.load(open(model,'rb'))
    # whether to perform normalization on results or not
    normalization = True if args.normalization=='Yes' else False
    # prefix of the output file
    prefix=args.prefix
    predictors = [x for x in my_data.columns if x!=id_column]
    # use the model to predict and test the performance
    activity=gbm_tuned.predict(my_data[predictors])
    if normalization:
        activity=perform_normalization(activity)
    result = pd.concat([my_data[id_column], pd.Series(activity, name='score')], axis=1)
    print ('')
    print (result.describe())
    print ('')
    result.to_csv('%s_result.txt'%(prefix), sep='\t', index=None)

if __name__ == '__main__':
    try:
        args=GBR_pred_parseargs()
        GBR_pred_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)
