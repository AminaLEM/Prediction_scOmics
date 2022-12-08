#!/usr/bin/env python
# coding: utf-8

# In[36]:


from sklearn.preprocessing import LabelEncoder, minmax_scale
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
from matplotlib.backends.backend_pdf import PdfPages
import sklearn.metrics as skm
import matplotlib.pyplot as plt
from collections import defaultdict
import shap
from numba import njit, prange
import seaborn as sns
from statannot import add_stat_annotation

def confidence_interval( values):
        mean = np.mean(values)
        alpha = 0.95
        p = ((1.0-alpha)/2.0) * 100
        bottom = max(0.0, np.percentile(values, p))
        p = (alpha+((1.0-alpha)/2.0)) * 100
        top = min(1.0, np.percentile(values, p))
        return mean, bottom, top
    
def eval_box(metrics,tit):
    f = plt.figure()
    val_vec=list()
    bt_vec=list()
    tp_vec=list()
    x_vec=list()
    for d in metrics:
        val, bt, tp= confidence_interval(metrics[d])
        val_vec.append(val)
        bt_vec.append(val-bt)
        tp_vec.append(tp-val)
        x_vec.append(d)
    
    with plt.rc_context({'figure.figsize': (4, 3), 'figure.dpi':300, "font.size" : 16}):
        plt.errorbar(x_vec, val_vec, yerr=(bt_vec, tp_vec),fmt='o', 
            capsize=5,
            ecolor='k', 
            lw=3,
            ls = ':',
            color='blue')
        plt.title(tit)
        plt.ylim(0.4, 1.05)
        plt.plot([], c='k', label='CI 95%')
        plt.plot([], c='blue', label='mean')
        plt.legend(loc="lower right")

    return f
path= snakemake.params["path_inp"]

def eval_box_sbn(metrics,tit):
    width=0.6
    sns.set(rc={"figure.figsize":(4, 5)})    
    f = plt.figure()
    data=pd.DataFrame([], columns=["modality",'value'])
    for d in metrics:
        data_=pd.DataFrame([np.full((len(metrics[d])), d),metrics[d]], index=["modality",'value'])
        data=pd.concat([data,data_.transpose()],ignore_index=True)
    data.value= data.value.astype('float')
    mod=data.modality.unique()
    if "GE" in mod:
        pairs=[("GE", "CC"),("CC", "CC&GE"),("CC&GE\n+GenderAge", "CC&GE")]

    else:
        pairs=[("CC&GE", "SVM"),("CC&GE", "RF"),("CC&GE", "LogisticR")]
    sns.set_style('whitegrid')
    with plt.rc_context({'figure.dpi':300, "font.size" : 16}):

        g = sns.boxplot(x="modality", y='value', data=data, width=width)

        add_stat_annotation(g, x="modality", y='value', data=data,
                        box_pairs=pairs,
                        test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
        plt.ylabel(tit)
        if len(mod)==3:
            plt.xlabel('modality')
        else:
            plt.xlabel('model')
    return f


def compute_metrics(y):
    metrics = defaultdict(list)
    col_len= y.shape[1]
    for i in range(col_len):
            l=y.columns[i]
            metrics['auc'].append(skm.roc_auc_score(Ytest,y[l]))
            metrics['acc'].append(skm.accuracy_score(Ytest,y[l]>=0.5))
            metrics['f1'].append(skm.f1_score(Ytest,y[l]>=0.5))
            metrics['rec'].append(skm.recall_score(Ytest,y[l]>=0.5))
            metrics['prc'].append(skm.precision_score(Ytest,y[l]>=0.5))
            metrics['auprc'].append(skm.average_precision_score(Ytest,y[l]))
    return metrics

def predict_loop(model_,data):
    y_score1 = pd.DataFrame([])    
    model_len=len(model_)
    for i in range(model_len):
        model=model_[i]
        y_score1['sampling'+str(i)]=model.predict(data).flatten()
    return y_score1

# path to output of snakemake model (pred folder + extension see below)
out_fig = path+'/pred/MLP'
out_fig_GA = path.replace('output_Top5','output_Top5_GenderAge')+'/pred/MLP'
svm_fig=path+'/pred/SVM'
LogReg_fig=path+'/pred/LogisticR'
RF_fig=path+'/pred/RF'


# # PLot meterics compared to baselines and Gender age CC, GE and CC&GE

# In[37]:


# this is just to load ground truth Ytest
x_cell=pd.DataFrame([])
for elem in [path+'/korean.h5Seurat/annotation.csv'
             , path+'/stanford.h5Seurat/annotation.csv']:
    x_=pd.read_csv (elem,index_col=0)
    x_cell = pd.concat([x_cell, x_], axis =0)

label= x_cell.iloc[:,-1].values
le = LabelEncoder()
Ytest = le.fit_transform(label)


# In[38]:


MLP_j= compute_metrics(pd.read_csv(out_fig + '_CC_GE.csv'))
MLP_CC= compute_metrics(pd.read_csv(out_fig + '_CC.csv'))
MLP_GE= compute_metrics(pd.read_csv(out_fig + '_GE.csv'))
MLP_j_GA= compute_metrics(pd.read_csv(out_fig_GA + '_CC_GE.csv'))

SVM_j= compute_metrics(pd.read_csv(svm_fig + '_CC_GE.csv'))
SVM_CC= compute_metrics(pd.read_csv(svm_fig + '_CC.csv'))
SVM_GE= compute_metrics(pd.read_csv(svm_fig + '_GE.csv'))

LogReg_j= compute_metrics(pd.read_csv(LogReg_fig + '_CC_GE.csv'))
LogReg_CC= compute_metrics(pd.read_csv(LogReg_fig + '_CC.csv'))
LogReg_GE= compute_metrics(pd.read_csv(LogReg_fig + '_GE.csv'))

RF_j= compute_metrics(pd.read_csv(RF_fig + '_CC_GE.csv'))
RF_CC= compute_metrics(pd.read_csv(RF_fig + '_CC.csv'))
RF_GE= compute_metrics(pd.read_csv(RF_fig + '_GE.csv'))


# In[40]:


#Comparison to baseline 
fig14=eval_box_sbn({'LogisticR':LogReg_j['auc'],'SVM':SVM_j['auc'], 'RF':RF_j['auc'], 'CC&GE':MLP_j['auc']},'AUC')
fig15=eval_box_sbn({'LogisticR':LogReg_j['auprc'],'SVM':SVM_j['auprc'], 'RF':RF_j['auprc'], 'CC&GE':MLP_j['auprc']},'Average Precision')
fig16=eval_box_sbn({'LogisticR':LogReg_j['acc'],'SVM':SVM_j['acc'], 'RF':RF_j['acc'], 'CC&GE':MLP_j['acc']},'Accuracy')
fig17=eval_box_sbn({'LogisticR':LogReg_j['prc'],'SVM':SVM_j['prc'], 'RF':RF_j['prc'], 'CC&GE':MLP_j['prc']},'Precision')
fig18=eval_box_sbn({'LogisticR':LogReg_j['rec'],'SVM':SVM_j['rec'], 'RF':RF_j['rec'], 'CC&GE':MLP_j['rec']},'Recall')
fig19=eval_box_sbn({'LogisticR':LogReg_j['f1'],'SVM':SVM_j['f1'], 'RF':RF_j['f1'], 'CC&GE':MLP_j['f1']},'F1-Score')
pp = PdfPages(path+'/figures/baseline.pdf')
pp.savefig(fig14, bbox_inches='tight')
pp.savefig(fig15, bbox_inches='tight')
pp.savefig(fig16, bbox_inches='tight')
pp.savefig(fig17, bbox_inches='tight')
pp.savefig(fig18, bbox_inches='tight')
pp.savefig(fig19, bbox_inches='tight')
pp.close()


# In[41]:


# Comaprison to CC, GE and GenderAge
fig8=eval_box_sbn({'GE':MLP_GE['auc'],'CC':MLP_CC['auc'],'CC&GE':MLP_j['auc'],'CC&GE\n+GenderAge':MLP_j_GA['auc']},'AUC')
fig9=eval_box_sbn({'GE':MLP_GE['auprc'],'CC':MLP_CC['auprc'],'CC&GE':MLP_j['auprc'],'CC&GE\n+GenderAge':MLP_j_GA['auprc']},'Average Precision')
fig10=eval_box_sbn({'GE':MLP_GE['acc'],'CC':MLP_CC['acc'],'CC&GE':MLP_j['acc'],'CC&GE\n+GenderAge':MLP_j_GA['acc']},'Accuracy')
fig11=eval_box_sbn({'GE':MLP_GE['prc'],'CC':MLP_CC['prc'],'CC&GE':MLP_j['prc'],'CC&GE\n+GenderAge':MLP_j_GA['prc']},'Precision')
fig12=eval_box_sbn({'GE':MLP_GE['rec'],'CC':MLP_CC['rec'],'CC&GE':MLP_j['rec'],'CC&GE\n+GenderAge':MLP_j_GA['rec']},'Recall')
fig13=eval_box_sbn({'GE':MLP_GE['f1'],'CC':MLP_CC['f1'],'CC&GE':MLP_j['f1'],'CC&GE\n+GenderAge':MLP_j_GA['f1']},'F1-Score')
pp = PdfPages(path+'/figures/predGECC.pdf')
pp.savefig(fig8, bbox_inches='tight')
pp.savefig(fig9, bbox_inches='tight')
pp.savefig(fig10, bbox_inches='tight')
pp.savefig(fig11, bbox_inches='tight')
pp.savefig(fig12, bbox_inches='tight')
pp.savefig(fig13, bbox_inches='tight')
pp.close()


# # PLot AUC per DEG pseudo_bulk

# In[11]:


from sklearn.preprocessing import LabelEncoder, minmax_scale
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
from matplotlib.backends.backend_pdf import PdfPages
import sklearn.metrics as skm
import matplotlib.pyplot as plt
from collections import defaultdict
import shap
from numba import njit, prange
import seaborn as sns
from statannot import add_stat_annotation

def confidence_interval( values):
        mean = np.mean(values)
        alpha = 0.95
        p = ((1.0-alpha)/2.0) * 100
        bottom = max(0.0, np.percentile(values, p))
        p = (alpha+((1.0-alpha)/2.0)) * 100
        top = min(1.0, np.percentile(values, p))
        return mean, bottom, top
    
def eval_box(metrics,tit):
    f = plt.figure()
    val_vec=list()
    bt_vec=list()
    tp_vec=list()
    x_vec=list()
    for d in metrics:
        val, bt, tp= confidence_interval(metrics[d])
        val_vec.append(val)
        bt_vec.append(val-bt)
        tp_vec.append(tp-val)
        x_vec.append(d)
    
    with plt.rc_context({'figure.figsize': (4, 3), 'figure.dpi':300, "font.size" : 16}):
        plt.errorbar(x_vec, val_vec, yerr=(bt_vec, tp_vec),fmt='o', 
            capsize=5,
            ecolor='k', 
            lw=3,
            ls = ':',
            color='blue')
        plt.title(tit)
        plt.ylim(0.4, 1.05)
        plt.plot([], c='k', label='CI 95%')
        plt.plot([], c='blue', label='mean')
        plt.legend(loc="lower right")

    return f

def eval_box_sbn(metrics,tit):
    f = plt.figure()
    data=pd.DataFrame([], columns=["modality",'value'])
    for d in metrics:
        data_=pd.DataFrame([np.full((len(metrics[d])), d),metrics[d]], index=["modality",'value'])
        data=pd.concat([data,data_.transpose()],ignore_index=True)
    data.value= data.value.astype('float')
    mod=data.modality.unique()
    if len(mod)==5:
        pairs=[("GE", "CC"),("GE", "CC&GE"),("CC", "CC&GE"),("CC&GE_GA", "CC&GE"),("CC&GE_10", "CC&GE")]
    else:
        pairs=[("ALL", "Top15"),("Top15", "Top10"),("Top5", "Top10")]
    sns.set(rc={"figure.figsize":(4, 5)})
    sns.set_style('whitegrid')
    with plt.rc_context({'figure.dpi':96, "font.size" : 16}):

        g = sns.boxplot(x="modality", y='value', data=data, width=0.6)
        
        add_stat_annotation(g, x="modality", y='value', data=data,
                        box_pairs=pairs,
                        test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
        plt.ylabel(tit)
        if len(mod)==3:
            plt.xlabel('modality')
        else:
            plt.xlabel('model')
    
    return f


def compute_metrics(y):
    metrics = defaultdict(list)
    col_len= y.shape[1]
    for i in range(col_len):
            l=y.columns[i]
            metrics['auc'].append(skm.roc_auc_score(Ytest,y[l]))
            metrics['acc'].append(skm.accuracy_score(Ytest,y[l]>=0.5))
            metrics['f1'].append(skm.f1_score(Ytest,y[l]>=0.5))
            metrics['rec'].append(skm.recall_score(Ytest,y[l]>=0.5))
            metrics['prc'].append(skm.precision_score(Ytest,y[l]>=0.5))
            metrics['auprc'].append(skm.average_precision_score(Ytest,y[l]))
    return metrics

def predict_loop(model_,data):
    y_score1 = pd.DataFrame([])    
    model_len=len(model_)
    for i in range(model_len):
        model=model_[i]
        y_score1['sampling'+str(i)]=model.predict(data).flatten()
    return y_score1


top5 = path+'/pred/MLP'
top10 = path.replace('output_Top5','output_Top10')+'/pred/MLP'
top15 = path.replace('output_Top5','output_Top15')+'/pred/MLP'
all_ = path.replace('output_Top5','output_All')+'/pred/MLP'

all__j= compute_metrics(pd.read_csv(all_ + '_CC_GE.csv'))
all__CC= compute_metrics(pd.read_csv(all_ + '_CC.csv'))
all__GE= compute_metrics(pd.read_csv(all_ + '_GE.csv'))

top10_j= compute_metrics(pd.read_csv(top10 + '_CC_GE.csv'))
top10_CC= compute_metrics(pd.read_csv(top10 + '_CC.csv'))
top10_GE= compute_metrics(pd.read_csv(top10 + '_GE.csv'))

top5_j= compute_metrics(pd.read_csv(top5 + '_CC_GE.csv'))
top5_CC= compute_metrics(pd.read_csv(top5 + '_CC.csv'))
top5_GE= compute_metrics(pd.read_csv(top5 + '_GE.csv'))

top15_j= compute_metrics(pd.read_csv(top15 + '_CC_GE.csv'))
top15_CC= compute_metrics(pd.read_csv(top15 + '_CC.csv'))
top15_GE= compute_metrics(pd.read_csv(top15 + '_GE.csv'))


# In[13]:


fig14=eval_box_sbn({'Top5':top5_GE['auc'],'Top10':top10_GE['auc'], 'Top15':top15_GE['auc'], 'ALL':all__GE['auc']},'AUC')
fig15=eval_box_sbn({'Top5':top5_GE['auprc'],'Top10':top10_GE['auprc'], 'Top15':top15_GE['auprc'], 'ALL':all__GE['auprc']},'Average Precision')
fig16=eval_box_sbn({'Top5':top5_GE['acc'],'Top10':top10_GE['acc'], 'Top15':top15_GE['acc'], 'ALL':all__GE['acc']},'Accuracy')
fig17=eval_box_sbn({'Top5':top5_GE['prc'],'Top10':top10_GE['prc'], 'Top15':top15_GE['prc'], 'ALL':all__GE['prc']},'Precision')
fig18=eval_box_sbn({'Top5':top5_GE['rec'],'Top10':top10_GE['rec'], 'Top15':top15_GE['rec'], 'ALL':all__GE['rec']},'Recall')
fig19=eval_box_sbn({'Top5':top5_GE['f1'],'Top10':top10_GE['f1'], 'Top15':top15_GE['f1'], 'ALL':all__GE['f1']},'F1-Score')
pp = PdfPages(path+'/figures/GEpred.pdf')
pp.savefig(fig14, bbox_inches='tight')
pp.savefig(fig15, bbox_inches='tight')
pp.savefig(fig16, bbox_inches='tight')
pp.savefig(fig17, bbox_inches='tight')
pp.savefig(fig18, bbox_inches='tight')
pp.savefig(fig19, bbox_inches='tight')
pp.close()


# In[ ]:





# In[ ]:




