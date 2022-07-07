# import scipy
from sklearn.preprocessing import LabelEncoder, minmax_scale
import pickle
import numpy as np
import pandas as pd
import random
import os
####*IMPORANT*: Have to do this line *before* importing tensorflow
import tensorflow as tf
from keras.callbacks import ModelCheckpoint
from sklearn.utils import resample
from keras.models import load_model
# from keras import backend as K
# print(tf.config.list_physical_devices("GPU"))

def reset_random_seeds():
   os.environ['PYTHONHASHSEED']=str(1)
   tf.random.set_seed(0)
   np.random.seed(1234)


x_cell=pd.read_csv (snakemake.input[0],index_col=0)
x_exp=pd.read_csv (snakemake.input[1],index_col=0)
x_exp=x_exp.loc[x_exp['condition'].isin(['Mild','Severe']),:]
x_cell=x_cell.loc[x_cell['condition'].isin(['Mild','Severe']),:]
x_cell=x_cell.drop(['Doublet','Eryth','NK_CD56bright'],axis=1)

train_set_f=snakemake.output[0]
val_set_f=snakemake.output[1]
model_j_f=snakemake.output[2]
model_e_f=snakemake.output[3]
model_c_f=snakemake.output[4]

path=snakemake.params[0]
x_cell=x_cell.loc[x_exp.index,:]

label= x_cell.iloc[:,-1].values
who=x_exp.iloc[:,-1].values
x_cell= x_cell.drop('condition',axis=1)
x_exp= x_exp.drop('condition',axis=1)
x_exp= x_exp.drop('who_score',axis=1)

genes = x_exp.columns


le = LabelEncoder()
Ytrain = le.fit_transform(label)

x_exp = minmax_scale(x_exp, axis = 0)
x_cell= x_cell.div(x_cell.sum(axis=1), axis=0)


def training(model, trainDataOne,y, valid_set):
    reset_random_seeds()
    # Optimizer setting    
    optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)

    # Model compiling settings
    model.compile(optimizer=optimizer,
              loss=tf.keras.losses.BinaryCrossentropy(),
              metrics=['accuracy'])

    # A mechanism that stops training if the validation loss is not improving for more than n_idle_epochs.
    n_idle_epochs = 100
    earlyStopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=n_idle_epochs, min_delta=0.001)
    mc = ModelCheckpoint(model_j_f+'.h5', monitor='val_loss', mode='min', save_best_only=True)
    
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
    n_epochs = 2000
    history = model.fit(
      trainDataOne, y, batch_size=1,
      epochs=n_epochs, validation_data = valid_set, verbose=0, callbacks=[log_display,earlyStopping,mc])
    return model, history




def build_classifier2(inp1, inp2):
    reset_random_seeds()
    # define two sets of inputs
    inputA = tf.keras.layers.Input(shape=(inp1.shape[1],))
    inputB = tf.keras.layers.Input(shape=(inp2.shape[1],))
    # the first branch operates on the first input
    x = tf.keras.layers.Dropout(0.2, seed = 0)(inputA)
    x = tf.keras.layers.Dense(10, activation="relu",kernel_initializer= tf.keras.initializers.GlorotNormal(seed=1234))(x)
    x = tf.keras.layers.Dense(4, activation="relu", kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(x)
    x = tf.keras.Model(inputs=inputA, outputs=x)
    # the second branch opreates on the second input
    y = tf.keras.layers.Dropout(0.2, seed = 0)(inputB)
#     y = tf.keras.layers.Dense(10, activation="relu",
#                               kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(y)  
    y = tf.keras.layers.Dense(10, activation="relu", kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(y)  
    y = tf.keras.layers.Dense(4, activation="relu", kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(y)  

    #     y = tf.keras.layers.BatchNormalization()(y)        
    y = tf.keras.Model(inputs=inputB, outputs=y)
    # combine the output of the two branches
    combined = tf.keras.layers.concatenate([x.output, y.output])
#     z = tf.keras.layers.Dense(4, activation="relu"                              
#                               , kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(combined)
# #     z = tf.keras.layers.Dense(10, activation="relu"                              
#                               , kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(combined)
    z = tf.keras.layers.Dense(2, activation="relu"                              
                              , kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(combined)
   
    z = tf.keras.layers.Dense(1, activation="sigmoid", 
                              kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234))(z)
    # our model will accept the inputs of the two branches and
    # then output a single value
    model = tf.keras.Model(inputs=[x.input, y.input], outputs=z)
    return model

# In[42]:

# In[42]:


def build_classifier(inp1):
    model = tf.keras.Sequential([
    tf.keras.layers.Dropout(0.2, seed = 0,input_shape=(inp1.shape[1],)),
    tf.keras.layers.Dense(10, activation='relu', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234)),  
    tf.keras.layers.Dense(4, activation='relu', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234)),        
    tf.keras.layers.Dense(2, activation='relu', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234)),        
        
    tf.keras.layers.Dense(1, activation='sigmoid', kernel_initializer=tf.keras.initializers.GlorotNormal(seed=1234)),
    ])
    return model


# In[43]:


def provide_stratified_bootstap_sample_indices(bs_sample,percent):

    strata = pd.DataFrame(who).value_counts()
    bs_index_list_stratified = []

    for idx_stratum_var, n_stratum_var in strata.iteritems():
        data_index_stratum = list(np.where(who == idx_stratum_var[0])[0])
        kk=round(len(data_index_stratum )*percent)
        bs_index_list_stratified.extend(random.sample(data_index_stratum , k = kk))
    return bs_index_list_stratified


sets = np.concatenate((x_exp, x_cell), axis=1)
sets = np.column_stack((sets, Ytrain))
dim_exp = x_exp.shape[1]
dim_cells = x_cell.shape[1]
n_iterations = 30
train_set = list()
val_set = list()
model_j={}
model_e={}
model_c={}
random.seed(42)
seeds = random.sample(range(0, 500), 30)
for i in range(n_iterations):
    random.seed(seeds[i])
    print(i)
    # prepare train and test sets
    all_sets= x_exp,
    #train = resample(sets, n_samples=n_size)
    bs_index_list_stratified= provide_stratified_bootstap_sample_indices(sets,0.8)
    train= sets[bs_index_list_stratified , :]
    test = np.array([x for x in sets if x.tolist() not in train.tolist()])
    # fit model
    model = build_classifier2(train[:,:dim_exp],train[:,dim_exp:(dim_exp+dim_cells)])
    model, history= training(model, [train[:,:dim_exp],train[:,dim_exp:(dim_exp+dim_cells)]]
                             , train[:,-1],([test[:,:dim_exp],test[:,dim_exp:(dim_exp+dim_cells)]],test[:,-1]))
    model_j[i]=load_model(model_j_f+'.h5')
    # evaluate model
    model = build_classifier(train[:,:dim_exp])
    model, history= training(model, train[:,:dim_exp]
                             , train[:,-1],(test[:,:dim_exp],test[:,-1]))
    model_e[i]=load_model(model_j_f+'.h5')
 
    # evaluate model
    model = build_classifier(train[:,dim_exp:(dim_exp+dim_cells)])
    model, history= training(model, train[:,dim_exp:(dim_exp+dim_cells)]
                             , train[:,-1],(test[:,dim_exp:(dim_exp+dim_cells)],test[:,-1]))
    model_c[i]=load_model(model_j_f+'.h5')
    
    # evaluate model 
    val_set.append(test)    
    train_set.append(train)        


with open(val_set_f, 'wb') as b:
    pickle.dump(val_set,b)
with open(train_set_f, 'wb') as b:
    pickle.dump(train_set,b)
with open(model_j_f, 'wb') as b:
    pickle.dump(model_j,b)
with open(model_e_f, 'wb') as b:
    pickle.dump(model_e,b)
with open(model_c_f, 'wb') as b:
    pickle.dump(model_c,b)
    