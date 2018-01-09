import os
import sys
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp
import Bio.SeqUtils.MeltingTemp as Tm

#this function is used to convert the positon-independent feature that is represented by vector to the feature is represented by a value. 
def Vector_feature_to_Value_feature(position_feature_depedent_Dic,basepairLst):
    new_Dic={}
    for position in position_feature_depedent_Dic:
        for i,base in enumerate(basepairLst):
            new_Dic['%s_%s'%(position,base)]=position_feature_depedent_Dic[position][i]
    return new_Dic

#######################
#this function is used to extract the Order1 feature
def Order1(sequence):
    n=len(sequence)
    seq=sequence
    baseDic={}
    order=1
    baseLst=['A','T','C','G']
    position_independentDic={}
    position_dependentDic={}
#generate the Base and its vector.
    for i,base in enumerate(baseLst):
        baseDic[base]=np.zeros(4**order)
        baseDic[base][i]=1 
        position_independentDic['order1_IP_%s'%(base)]=0

#extract the seq feature
    for i in range(len(seq)):
        for j,base in enumerate(baseDic):
            if base==seq[i]:
                position_dependentDic['order1_P%s'%(i+1)]=baseDic[base]
                position_independentDic['order1_IP_%s'%(base)]+=1
#    position_dependent=sorted(position_dependent.items(),key=lambda item:item[0])
#    position_independent=sorted(position_independent.items(),key=lambda item:item[0])
    position_dependentDic=Vector_feature_to_Value_feature(position_dependentDic,baseLst)
    Order1_positionDic=dict(position_dependentDic.items()+position_independentDic.items())
    return Order1_positionDic
############################
#this fuction is used to extract the Order2 feature
def Order2(sequence):
    seq=sequence
    BasepairDic={}
    BasepairLst=[]
    position_dependentDic={}
    position_independentDic={}
    order=2
    baseLst=['A','T','C','G']
#generate the Basepair and its vector that contain 0 or 1.
    for base1 in baseLst:
        for base2 in baseLst:
            BasepairLst.append(base1+base2)
    for i in range(len(BasepairLst)):
        BasepairDic[BasepairLst[i]]=np.zeros(4**order)
        BasepairDic[BasepairLst[i]][i]=1
        position_independentDic['order2_IP_%s'%(BasepairLst[i])]=0
#extract the seqence feature
    for i in range(len(seq)-1):
        seq_pair=seq[i:i+2]
        for j,basepair in enumerate(BasepairLst):
            if seq_pair==basepair:
                position_dependentDic['order2_P%s'%(i+1)]=BasepairDic[basepair]
                position_independentDic['order2_IP_%s'%(basepair)]+=1
#    position_dependent=sorted(position_dependent.items(),key=lambda item:item[0])
#    position_independent=sorted(position_independent.items(),key=lambda item:item[0])
#    print(Basepair)
    position_dependentDic=Vector_feature_to_Value_feature(position_dependentDic,BasepairLst)
    Order2_positonDic=dict(position_dependentDic.items()+position_independentDic.items())
    return Order2_positonDic

############################
#this function is used to extract the NGGN feture
def NGGNfeature(sequence):
    NNDic={}
    seq=sequence
    BasepairDic={}
    BasepairLst=[]
    baseLst=['A','T','G','C']
    for base1 in baseLst:
        for base2 in baseLst:
            BasepairLst.append(base1+base2)
    for i in range(len(BasepairLst)):
        BasepairDic[BasepairLst[i]]=np.zeros(4**2)
        BasepairDic[BasepairLst[i]][i]=1
    for basepair in BasepairLst:
        if basepair == seq:
            NNDic['NGGN']=BasepairDic[basepair]
    NNDic=Vector_feature_to_Value_feature(NNDic,BasepairLst)
    return NNDic

####################
#this function is used to extract the melt tempreture of sequence
def Temper(sequence):
    seq=sequence
    seq_7=seq[:7]
    seq_8=seq[7:15]
    seq_5=seq[15:20]
    TDic={}
    TDic['T20']=Tm.Tm_staluc(seq)
    TDic['T7']=Tm.Tm_staluc(seq_7)
    TDic['T8']=Tm.Tm_staluc(seq_8)
    TDic['T5']=Tm.Tm_staluc(seq_5)
    return TDic

##############################    
#this function is used to get the Reverse complementary sequence!
def complement(sequence):
    seq=''
    for rec in sequence:
        if rec=='A':
            rec='T'
        elif rec=='T':
            rec='A'
        elif rec=='C':
            rec='G'
        else:
            rec='C'
        seq+=rec
    seq=seq[::-1]
    return seq

########################
#this fuction is used to extract the all features of sgRNA sequence!
def feature(sequence,NGGN):
    Order1Position=Order1(sequence)
    Order2Position=Order2(sequence)
    Temprature=Temper(sequence)
    NGGN_sequence=NGGNfeature(NGGN)
    seq_feature=dict(Order1Position.items()+Order2Position.items()+Temprature.items()+NGGN_sequence.items())
    seq_feature['GC']=float(GC(sequence))/100
    seq_feature_name=sorted(Order1Position.keys())+sorted(Order2Position.keys())+sorted(Temprature.keys())+sorted(NGGN_sequence.keys())
    seq_feature_name.append('GC')
    return seq_feature_name,seq_feature

#######################################################################


######################################################################
def get_feature_main():
    sgRNA_fasta_file=sys.argv[1]
    sgRNAfasta=SeqIO.parse(sgRNA_fasta_file,"fasta")
    prefix=sys.argv[2]
    sgRNAfastaDic={}
    sgrnaFeatureDic={}
    ###########################################
    # check every sgRNA in the fasta file
    for rec in sgRNAfasta:
        sgRNAfastaDic[rec.id]=str(rec.seq)
        sgrnaFeatureDic[rec.id]={}
    ###########################################
    featureName=[]
    for sgrna_name in sgRNAfastaDic:
        sgrnaSequence=sgRNAfastaDic[sgrna_name][4:24]
        NGGNSequence=sgRNAfastaDic[sgrna_name][24]+sgRNAfastaDic[sgrna_name][27]
        featureName,featureValue=feature(sgrnaSequence,NGGNSequence)
        sgrnaFeatureDic[sgrna_name]=featureValue
    output='%s_feature.txt'%(prefix)
    os.system('cat /dev/null >%s'%(output))
    with open(output,'w') as f1:
        writtenLine='sgRNAID,'
        for featureID in featureName:
            writtenLine+='%s,'%(featureID)
        f1.write(writtenLine[:-1]+'\n')
        for sgrna_name in sgRNAfastaDic:
            writtenLine='%s,'%(sgrna_name)
            for featureID in featureName:
                writtenLine+='%s,'%(sgrnaFeatureDic[sgrna_name][featureID])
            f1.write(writtenLine[:-1]+'\n')

if __name__ =='__main__':
    try:
        get_feature_main()
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)
