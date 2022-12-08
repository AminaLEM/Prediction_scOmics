#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import sklearn.metrics as skm
from sklearn.preprocessing import LabelEncoder, minmax_scale
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


# # PLot AUC per DEG of every celltype

# In[5]:


celltypes = ['gd T cells', 'Intermediate monocytes','Plasmablasts']
data_test= {}
data_train= {}
y_test= {}
y_train= {}
path= snakemake.params["path_inp"]
path_celltype= snakemake.params["path_inp_celltype"]

for ct in celltypes:
    data_1= pd.read_csv (path_celltype+'/stanford.h5Seurat/'+ct+'/selected_ge.csv',index_col=0)
    data_2= pd.read_csv (path_celltype+'/korean.h5Seurat/'+ct+'/selected_ge.csv',index_col=0)
    ct_=ct+" ("+data_1.columns[0]+")"
    data_test[ct_] = pd.concat([data_1, data_2], axis =0)  
    data_train[ct_]=pd.read_csv (path_celltype+'/merged_training/'+ct+'/selected_ge.csv',index_col=0)
    Y_test= data_test[ct_].condition
    Y_train= data_train[ct_].condition
    data_test[ct_] = data_test[ct_][data_test[ct_].columns[0]]/ data_test[ct_][data_test[ct_].columns[0]].max()
    data_train[ct_] = data_train[ct_][data_train[ct_].columns[0]]/ data_train[ct_][data_train[ct_].columns[0]].max() 
    le = LabelEncoder()
    y_test[ct_] = le.fit_transform(Y_test)    
    y_train[ct_] = le.fit_transform(Y_train)    


# In[6]:


with plt.rc_context({'figure.figsize': (6, 5), 'figure.dpi':300, "font.size" : 10}):
            for key in data_test.keys():
                ns_auc = skm.roc_auc_score(y_test[key], data_test[key])
                # calculate roc curves
                lr_fpr, lr_tpr, _ = skm.roc_curve(y_test[key], data_test[key])
                # plot the roc curve for the model
                plt.plot(lr_fpr, lr_tpr, marker='.', label=key+": "+str(round(ns_auc,2)))
                # axis labels
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                # show the legend
                plt.legend()
                # show the plot
plt.show()

with plt.rc_context({'figure.figsize': (6, 5), 'figure.dpi':300, "font.size" : 10}):
            for key in data_test.keys():

                lr_precision, lr_recall, _ = skm.precision_recall_curve(y_test[key], data_test[key])
                lr_auc = skm.auc(lr_recall, lr_precision)
                avr_prec= skm.average_precision_score(y_test[key], data_test[key])
                # plot the precision-recall curves
                plt.plot(lr_recall, lr_precision, marker='.', label=key+": "+str(round(lr_auc,2))+"_"+str(round(avr_prec,2)))
                # axis labels
                plt.xlabel('Recall')
                plt.ylabel('Precision')
                # show the legend
                plt.legend()
plt.show()


# In[7]:


AUC = pd.DataFrame([])
for key in data_test.keys():
                AUC = pd.concat([AUC, pd.DataFrame([key,'AUC',skm.roc_auc_score(y_test[key], data_test[key])]).transpose()], axis =0)

for key in data_test.keys():
                AUC = pd.concat([AUC, pd.DataFrame([key,'Average Precision',skm.average_precision_score(y_test[key], data_test[key])]).transpose()], axis =0)
AUC = pd.concat([AUC, pd.DataFrame(['Classical monocytes','AUC',0.73,]).transpose()], axis =0)   
AUC = pd.concat([AUC, pd.DataFrame(['Classical monocytes','Average Precision',0.83,]).transpose()], axis =0)               

AUC.columns = ['DEG per celltype','Meteric','Value']
AUC_ = AUC.pivot(index = 'DEG per celltype',columns='Meteric', values='Value')
f = plt.figure()
with plt.rc_context({'figure.figsize': (6, 5), 'figure.dpi':300, "font.size" : 10}):
    AUC_.plot(kind='barh',color=['#3693a4', '#f7464e'])
    plt.legend(fontsize=7, loc='center right')
    plt.savefig(path+'/figures/AUC.PR_DEcelltypes.pdf',dpi=300, bbox_inches='tight')


# # PLot AUC per celltype proportion

# In[8]:


stan= pd.read_csv (path+'/stanford.h5Seurat/annotation.csv',index_col=0)
kor= pd.read_csv (path+'/korean.h5Seurat/annotation.csv',index_col=0)
bonn_berlin= pd.read_csv (path+'/merged_training/annotation.csv',index_col=0)
test=pd.concat([stan, kor], axis =0)

bonn_berlin.columns.values[bonn_berlin.columns == 'Progenitor cells']= "Platelet cells"
test.columns.values[test.columns == 'Progenitor cells']= "Platelet cells"


# In[9]:


data_test= {}
data_train= {}
y_test= {}
y_train= {}
celltypes= ['Classical monocytes',"Non classical monocytes","Intermediate monocytes", 
            "Plasmacytoid dendritic cells","Myeloid dendritic cells","Plasmablasts", "B cells",
             "CD8 T cells" , "CD4 T cells","T regulatory cells", "MAIT cells", "gd T cells",
            "Natural killer cells","Low-density neutrophils" ,"Platelet cells"]
#celltypes = stan.columns[stan.columns!='condition']
for ct in celltypes:
    data_test[ct] = test[ct]  
    data_train[ct]=bonn_berlin[ct]
    data_test[ct] = data_test[ct]/ test.drop('condition', axis=1).sum(axis=1)
    data_train[ct] = data_train[ct]/ bonn_berlin.drop('condition', axis=1).sum(axis=1)
le = LabelEncoder()
y_test= le.fit_transform(test.condition)    
y_train = le.fit_transform(bonn_berlin.condition)  
    


# In[10]:


with plt.rc_context({'figure.figsize': (6, 5), 'figure.dpi':300, "font.size" : 10}):
            for key in data_test.keys():
                ns_auc = skm.roc_auc_score(y_test, data_test[key])
                # calculate roc curves
                lr_fpr, lr_tpr, _ = skm.roc_curve(y_test, data_test[key])
                # plot the roc curve for the model
                plt.plot(lr_fpr, lr_tpr, marker='.', label=key+": "+str(round(ns_auc,2)))
                # axis labels
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                # show the legend
                plt.legend()
                # show the plot
plt.show()

with plt.rc_context({'figure.figsize': (6, 5), 'figure.dpi':300, "font.size" : 10}):
            for key in data_test.keys():

                lr_precision, lr_recall, _ = skm.precision_recall_curve(y_test, data_test[key])
                lr_auc = skm.auc(lr_recall, lr_precision)
                avr_prec= skm.average_precision_score(y_test, data_test[key])
                # plot the precision-recall curves
                plt.plot(lr_recall, lr_precision, marker='.', label=key+": "+str(round(lr_auc,2))+"_"+str(round(avr_prec,2)))
                # axis labels
                plt.xlabel('Recall')
                plt.ylabel('Precision')
                # show the legend
                plt.legend()
plt.show()


# In[11]:


AUC = pd.DataFrame([])
for key in data_test.keys():
                AUC = pd.concat([AUC, pd.DataFrame([key,'AUC',skm.roc_auc_score(y_test, data_test[key])]).transpose()], axis =0)

for key in data_test.keys():
                AUC = pd.concat([AUC, pd.DataFrame([key,'Average Precision',skm.average_precision_score(y_test, data_test[key])]).transpose()], axis =0)
                
AUC.columns = ['Cell types','Meteric','Value']
AUC_ = AUC.pivot(index = 'Cell types',columns='Meteric', values='Value')
AUC_= AUC_.loc[celltypes,:]
f = plt.figure()
with plt.rc_context({'figure.figsize': (6, 5), 'figure.dpi':300, "font.size" : 10}):
    AUC_.plot(kind='barh',color=['#3693a4', '#f7464e'])
    plt.legend(fontsize=7, loc='center right')
    plt.savefig(path+'/figures/AUC.PR_celltypes.pdf',dpi=300, bbox_inches='tight')


# In[ ]:




