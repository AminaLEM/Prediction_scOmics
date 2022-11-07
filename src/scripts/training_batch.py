# import scipy
from sklearn.preprocessing import LabelEncoder, minmax_scale
import pickle
import numpy as np
import pandas as pd
import random
import os
import tensorflow as tf
from keras.callbacks import ModelCheckpoint
from sklearn.utils import resample
from keras.models import load_model
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2,f_classif

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from collections import defaultdict
from hyperopt import  SparkTrials, STATUS_OK, tpe,fmin, hp
from joblib import Parallel, delayed
# from numba import njit, prange
# from keras import backend as K
# print(tf.config.list_physical_devices("GPU"))
x_cell=pd.read_csv (snakemake.input['CC'],index_col=0)
x_exp=pd.read_csv (snakemake.input['GE'],index_col=0)
train_set_f=snakemake.output['train_set']
val_set_f=snakemake.output['val_set']
model_j_f=snakemake.output['model_j']
model_e_f=snakemake.output['model_e']
model_c_f=snakemake.output['model_c']
svm_e_f=snakemake.output['svm_e']
svm_c_f=snakemake.output['svm_c']
LReg_e_f=snakemake.output['LReg_e']
LReg_c_f=snakemake.output['LReg_c']
LogReg_e_f=snakemake.output['LogReg_e']
LogReg_c_f=snakemake.output['LogReg_c']
# Lasso_e_f=snakemake.output[11]
# Lasso_c_f=snakemake.output[12]
RF_e_f=snakemake.output['RF_e']
RF_c_f=snakemake.output['RF_c']
svm_j_f=snakemake.output['svm_j']
LReg_j_f=snakemake.output['LReg_j']
LogReg_j_f=snakemake.output['LogReg_j']
# Lasso_j_f=snakemake.output[18]
RF_j_f=snakemake.output['RF_j']
path=snakemake.params[0]
from tensorflow.keras import regularizers
def reset_random_seeds():
   os.environ['PYTHONHASHSEED']=str(1)
   tf.random.set_seed(0)
   np.random.seed(1234)


def training(model, trainDataOne,y, valid_set, ref,space_):
    reset_random_seeds()
    # Optimizer setting    
    optimizer = tf.keras.optimizers.Adam(learning_rate=space_['lr'])

    # Model compiling settings
    model.compile(optimizer=optimizer,
              loss=tf.keras.losses.BinaryCrossentropy(),
              metrics=['accuracy'])

    # A mechanism that stops training if the validation loss is not improving for more than n_idle_epochs.
    n_idle_epochs = 100
    earlyStopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=n_idle_epochs, min_delta=0.001)
    mc = ModelCheckpoint(model_j_f+ref+'.h5', monitor='val_loss', mode='min', save_best_only=True)
    
    # Creating a custom callback to print the log after a certain number of epochs
    class NEPOCHLogger(tf.keras.callbacks.Callback):
        def __init__(self,per_epoch=100):
            '''
            display: Number of batches to wait before outputting loss
            '''
            self.seen = 0
            self.per_epoch = per_epoch

        def on_epoch_end(self, epoch, logs=None):
           if epoch % self.per_epoch == 0:
            print('Epoch {}, loss {:.2f}, val_loss {:.2f}, accuracy {:.2f}, val_accuracy {:.2f}'.format(epoch, logs['loss'], logs['val_loss'], logs['accuracy'], logs['val_accuracy']))

    log_display = NEPOCHLogger(per_epoch=100)
    # Training loop
    n_epochs = 3000
    history = model.fit(
      trainDataOne, y, batch_size=space_['batch'],
      epochs=n_epochs, validation_data = valid_set, verbose=0, callbacks=[log_display,earlyStopping,mc])
    return model, history




def build_j_classifier(inp1, inp2,space_):
   
    reset_random_seeds()
    # define two sets of inputs
    inputA = tf.keras.layers.Input(shape=(inp1.shape[1],))
    inputB = tf.keras.layers.Input(shape=(inp2.shape[1],))
    # the first branch operates on the first input
    x = tf.keras.layers.Dropout(space_['dropout'], seed = 0)(inputA)    
    x = tf.keras.layers.Dense(2, activation='relu',kernel_initializer= tf.keras.initializers.GlorotNormal(seed=1234))(x)
    x = tf.keras.Model(inputs=inputA, outputs=x)
    # the second branch opreates on the second input
    y = tf.keras.layers.Dropout(space_['dropout'], seed = 0)(inputB)
    y = tf.keras.layers.Dense(2, activation='relu', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(y)  
    y = tf.keras.Model(inputs=inputB, outputs=y)
    # combine the output of the two branches
    combined = tf.keras.layers.concatenate([x.output, y.output])
    z = tf.keras.layers.Dense(2, activation='relu', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(combined)   
    z = tf.keras.layers.Dense(1, activation="sigmoid", 
                              kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(z)
    # our model will accept the inputs of the two branches and
    # then output a single value
    model = tf.keras.Model(inputs=[x.input, y.input], outputs=z)
    return model


def build_classifier(inp1,space_):
    model = tf.keras.Sequential([
    tf.keras.layers.Dropout(space_['dropout'], seed = 0,input_shape=(inp1.shape[1],)),
        
    tf.keras.layers.Dense(2, activation='relu', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234)),        
        
    tf.keras.layers.Dense(1, activation='sigmoid', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234)),
    ])
    return model


def provide_stratified_bootstap_sample_indices(bs_sample,percent):

    strata = pd.DataFrame(who).value_counts()
    bs_index_list_stratified = []

    for idx_stratum_var, n_stratum_var in strata.iteritems():
        data_index_stratum = list(np.where(who == idx_stratum_var[0])[0])
        kk=round(len(data_index_stratum )*percent)
        bs_index_list_stratified.extend(random.sample(data_index_stratum , k = kk))
    return bs_index_list_stratified

# @njit(parallel=True)
def task(bs_list_stratified):
        bs_index_list_stratified, space_= bs_list_stratified
        res_task= {}
        train= sets[bs_index_list_stratified , :]
        res_task['train']= train
        test = np.array([x for x in sets if x.tolist() not in train.tolist()])
        res_task['test']=test
#         ref= str(id(res_task)+id(test))
        ref = str(random.random())
        # fit model
        y=train[:,-1]
        model = build_j_classifier(train[:,:dim_exp],train[:,dim_exp:(dim_exp+dim_cells)],space_)
        model, res_task['history']= training(model, [train[:,:dim_exp],train[:,dim_exp:(dim_exp+dim_cells)]]
                                 , y,([test[:,:dim_exp],test[:,dim_exp:(dim_exp+dim_cells)]],test[:,-1]),ref,space_)

        res_task['model_j']=load_model(model_j_f+ref+'.h5')
        train_set=train[:,:(dim_exp+dim_cells)]
        
        # evaluate model
        train_set=train[:,:dim_exp]
        model = build_classifier(train[:,:dim_exp],space_)
        model,_= training(model, train_set, y,(test[:,:dim_exp],test[:,-1]),ref,space_)
        res_task['model_e']=load_model(model_j_f+ref+'.h5')
                
        # evaluate model
        train_set=train[:,dim_exp:(dim_exp+dim_cells)]
        model = build_classifier(train[:,dim_exp:(dim_exp+dim_cells)],space_)
        model,_= training(model, train_set
                                 , y,(test[:,dim_exp:(dim_exp+dim_cells)],test[:,-1]),ref,space_)
        res_task['model_c']=load_model(model_j_f+ref+'.h5')
        print('finished task -------------------------------------------'+ref)
        return res_task

def train_loop(space):
        n_iterations = 30
        random.seed(1234)
        seeds = random.sample(range(0, 1000), n_iterations)
        stratified_all= list()
        for i in range(n_iterations):
            random.seed(seeds[i])
            stratified_all.append([provide_stratified_bootstap_sample_indices(sets,0.8),space])
        loss =0
        res_=defaultdict(list)
        results = Parallel(n_jobs=10)(delayed(task)(i) for i in stratified_all)
        for result in results:
                for key in result.keys():
                    res_[key].append(result[key])
                loss=min(result['history'].history['val_loss'])+loss
        return {'loss': loss, 'status': STATUS_OK, 'model':res_}

def task_comaparative(bs_list_stratified):
        bs_index_list_stratified= bs_list_stratified
        res_task= {}
        train= sets[bs_index_list_stratified , :]
        test = np.array([x for x in sets if x.tolist() not in train.tolist()])
        y=train[:,-1]
        train_set=train[:,:(dim_exp+dim_cells)]       
        models={'LReg_j':LinearRegression(),'LogReg_j':LogisticRegression(),'svm_j':svm.SVC(),'RF_j':RandomForestClassifier()}
        for key in models.keys():
            res_task[key] = models[key].fit(train_set,train[:,-1])
        # evaluate model
        train_set=train[:,:dim_exp]             
        models={'LReg_e':LinearRegression(),'LogReg_e':LogisticRegression(),'svm_e':svm.SVC(),'RF_e':RandomForestClassifier()}
        for key in models.keys():
            res_task[key] = models[key].fit(train_set,y)                
        # evaluate model
        train_set=train[:,dim_exp:(dim_exp+dim_cells)]
        models={'LReg_c':LinearRegression(),'LogReg_c':LogisticRegression(),'svm_c':svm.SVC(),'RF_c':RandomForestClassifier()}
        for key in models.keys():
            res_task[key] = models[key].fit(train_set,y)
        return res_task
    
def compare():
        n_iterations = 30
        random.seed(1234)
        seeds = random.sample(range(0, 1000), n_iterations)
        stratified_all= list()
        for i in range(n_iterations):
            random.seed(seeds[i])
            stratified_all.append(provide_stratified_bootstap_sample_indices(sets,0.8))

        loss =0
        res_=defaultdict(list)
        results = Parallel(n_jobs=10)(delayed(task_comaparative)(i) for i in stratified_all)
        
#         with Pool() as pool: 
            
#             for result in pool.imap(task, stratified_all):
                # evaluate model
        for result in results:
                for key in result.keys():
                    res_[key].append(result[key])
        return res_
    
if __name__ == "__main__":
    res=defaultdict(list)
    x_exp=x_exp.loc[x_exp['condition'].isin(['Mild','Severe']),:]
    x_cell=x_cell.loc[x_cell['condition'].isin(['Mild','Severe']),:]
    x_cell=x_cell.loc[x_exp.index,:]
    
    label= x_cell.iloc[:,-1].values
    who=x_exp.iloc[:,-1].values
    x_cell= x_cell.drop('condition',axis=1)
    x_exp= x_exp.drop('condition',axis=1)
    x_exp= x_exp.drop('who_score',axis=1)
    
    genes = x_exp.columns
    le = LabelEncoder()
    Ytrain = le.fit_transform(label)
    x_exp=x_exp/ 10
    x_cell= x_cell.div(x_cell.sum(axis=1), axis=0)  
#     selected_cols= variance_threshold(x_cell,0)
#     x_cell=x_cell.loc[:,selected_cols]
    selected_cols=SelectKBest(f_classif, k='all').fit(x_cell, label).get_feature_names_out()[SelectKBest(f_classif, k='all').fit(x_cell, label).pvalues_ <0.05]

    x_cell=x_cell.loc[:,selected_cols]
#     x_cell=x_cell[['Low-density neutrophils','Progenitor cells','Naive B cells','Myeloid dendritic cells','Plasmablasts','Non-switched memory B cells','Switched memory B cells','Natural killer cells','T regulatory cells','MAIT cells']]
    
    sets = np.concatenate((x_exp, x_cell), axis=1)
    sets = np.column_stack((sets, Ytrain))
    dim_exp = x_exp.shape[1]
    dim_cells = x_cell.shape[1]
    
    
#     with Pool() as pool:   
#         for result in pool.map(task, stratified_all):
#             # evaluate model
#             for key in result.keys():
#                 res[key].append(result[key])
    space_tf = {'lr': hp.choice('lr',[0.001,0.0001]),'dropout': hp.choice('dropout',[0,0.1,0.2]),'batch': hp.choice('batch',[1,16,32])
#                 ,'layers': hp.choice('layers',[{'layer':'two'}, {'layer': 2}]) 
            }
    trials =  SparkTrials(parallelism=3)
    best_classifier = fmin(train_loop, space_tf, algo=tpe.suggest, max_evals=12, trials=trials, rstate=np.random.default_rng(42))
    print(best_classifier)
    res = trials.results[np.argmin([r['loss'] for r in trials.results])]['model']
    
#     res=best_classifier.model
    
    pp = PdfPages(model_j_f+"history_loss.pdf")
    for his in res['history']:
        f = plt.figure()
        # summarize history for loss
        plt.plot(his.history['loss'])
        plt.plot(his.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='upper left')
        pp.savefig(f, bbox_inches='tight')
    pp.close()
    
    pp = PdfPages(model_j_f+"history_acc.pdf")
    for his in res['history']:
        f = plt.figure()
        # summarize history for loss
        plt.plot(his.history['accuracy'])
        plt.plot(his.history['val_accuracy'])
        plt.title('model accuracy')
        plt.ylabel('accuracy')
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='upper left')
        pp.savefig(f, bbox_inches='tight')
    pp.close()
    res_comparative=compare()
    res = {**res, **res_comparative}
    map_name_file = {'test':val_set_f,'train':train_set_f,'model_j':model_j_f ,'model_e':model_e_f,'model_c':model_c_f,'svm_j':svm_j_f,'svm_e':svm_e_f,'svm_c':svm_c_f,'RF_j':RF_j_f,'RF_e':RF_e_f,'RF_c':RF_c_f,'LReg_j':LReg_j_f,'LReg_e':LReg_e_f,'LReg_c':LReg_c_f,'LogReg_j':LogReg_j_f,'LogReg_e':LogReg_e_f,'LogReg_c':LogReg_c_f}
    for key in map_name_file.keys():
        with open(map_name_file[key], 'wb') as b:
            pickle.dump(res[key],b)       
#     with open(model_j_f+'trials', 'wb') as b:
#             pickle.dump(trials,b)           
    l_h5 = os.listdir(path)
    for item in l_h5:
        if item.endswith(".h5"):
            os.remove(os.path.join(path, item))            

    pickling_trials = dict()

    for k, v in trials.__dict__.items():
        if not k in ['_spark_context', '_spark']:
            pickling_trials[k] = v

    pickle.dump(pickling_trials, open(model_j_f+'pickling_trials.hyperopt', 'wb'))
    pickle.dump(selected_cols, open(model_j_f+'selected_cols.pkl', 'wb'))
            
#     with open(val_set_f, 'wb') as b:
#         pickle.dump(val_set_all,b)
#     with open(train_set_f, 'wb') as b:
#         pickle.dump(train_set_all,b)
#     with open(model_j_f, 'wb') as b:
#         pickle.dump(model_j_all,b)
#     with open(model_e_f, 'wb') as b:
#         pickle.dump(model_e_all,b)
#     with open(model_c_f, 'wb') as b:
#         pickle.dump(model_c_all,b)
#     with open(svm_c_f, 'wb') as b:
#         pickle.dump(svm_c_all,b)
#     with open(svm_e_f, 'wb') as b:
#         pickle.dump(svm_e_all,b)   
#     with open(LReg_e_f, 'wb') as b:
#         pickle.dump(LReg_e_all,b)        
#     with open(LReg_c_f, 'wb') as b:
#         pickle.dump(LReg_c_all,b)          
#     with open(LogReg_e_f, 'wb') as b:
#         pickle.dump(LogReg_e_all,b)        
#     with open(LogReg_c_f, 'wb') as b:
#         pickle.dump(LogReg_c_all,b)       
#     with open(RF_e_f, 'wb') as b:
#         pickle.dump(RF_e_all,b)        
#     with open(RF_c_f, 'wb') as b:
#         pickle.dump(RF_c_all,b)
        
#     with open(RF_j_f, 'wb') as b:
#         pickle.dump(RF_j_all,b)   
#     with open(LogReg_j_f, 'wb') as b:
#         pickle.dump(LogReg_j_all,b)       
#     with open(LReg_j_f, 'wb') as b:
#         pickle.dump(LReg_j_all,b)       
#     with open(svm_j_f, 'wb') as b:
#         pickle.dump(svm_j_all,b)       
