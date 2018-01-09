# this script is used to predict the activity of sgRNA given the fasta file of sgRNAs

import os
import sys
import matplotlib.pyplot as plt
import numpy as np 
import ConfigParser

# The function processing the configure file
def configprocess(text):
    config=ConfigParser.ConfigParser()
    config.read(text)
    if len(config.sections())>1:
        print('Make a mistake in the configure file') and os._exit(0) 
    section=''.join(config.sections())
    varlist={option:config.get(section,option) for option in config.options(section)}
    return varlist

# configure file check
def metricJudge(varlist):
    model    =True if varlist['model'] in ['Cas9', 'eSpCas9'] else False
    normalization    =True if varlist['normalization'] in ['Yes', 'No'] else False
    targetfasta  =True
    prefix       =True
    for item in varlist:
        if not eval(item) :
            print('The varaint %s in configure file is wrong, please input the correct one according to the user manual!!!!'%(item))
            print(varlist[item])
            os._exit(1)

# extract the feature
def extract_feature(raw_fasta, prefix):
    os.system('python sgRNA_features.py %s %s'%(raw_fasta, prefix))

# predict sgRNA activity
def predict(feature_file, id_column, model, normalization, prefix):
    os.system('python GBR_prediction.py -d %s -i %s -m %s -n %s -p %s'%(feature_file, id_column, model, normalization, prefix))

def main(configureFile):
    varlist= configprocess(configureFile)
    metricJudge(varlist)
    # featurization of sgRNA fasta file (N4N20NGGN3)
    extract_feature(varlist['targetfasta'], varlist['prefix'])
    feature_file=varlist['prefix']+'_feature.txt'
    # perform prediction
    id_column='sgRNAID'
    model='saved_model/Cas9_sgRNA_activity_GBR.pickle' if varlist['model']=='Cas9' else 'saved_model/eSpCas9_sgRNA_activity_GBR.pickle'
    predict(feature_file, id_column, model, varlist['normalization'], varlist['prefix'])
    print ('prediction process is finalized, please check your result: %s_result.txt'%(varlist['prefix']))

if __name__ =='__main__':
    main(sys.argv[1])
