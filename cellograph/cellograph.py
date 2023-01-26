import os
from tensorflow.keras.utils import to_categorical
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import shuffle
from sklearn.metrics import classification_report

from spektral.layers import GCNConv

from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dropout, Dense
from tensorflow.keras import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import TensorBoard, EarlyStopping
import tensorflow as tf
from tensorflow.keras.regularizers import l2

from collections import Counter
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

import phate
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from math import *
import meld
import anndata as ad
from anndata import AnnData
from scipy import stats
from scipy.stats import pearsonr
import scipy.stats
import pylab as P
import seaborn as sns
import scvelo as scv


def encode_label(labels):
    label_encoder = LabelEncoder()
    labels = label_encoder.fit_transform(labels)
    labels = to_categorical(labels)
    return labels, label_encoder.classes_

def limit_data(labels,limit=20,val_num=400,test_num=400):
    '''
    Get the index of train, validation, and test data
    '''
    label_counter = dict((l, 0) for l in labels)
    train_idx = []

    for i in range(len(labels)):
        label = labels[i]
        if label_counter[label]<limit:
            train_idx.append(i)
            label_counter[label]+=1
        
        if all(count == limit for count in label_counter.values()):
            break
    
    
    rest_idx = [x for x in range(len(labels)) if x not in train_idx]
    val_idx = rest_idx[:val_num]
    test_idx = rest_idx[val_num:(val_num+test_num)]
    return train_idx, val_idx,test_idx

def run_cellograph(
                   adata = None, 
                   labels = None,
                   neighbors = 15,
                   pcs = 60, 
                   epochs= 200, 
                   num_classes = None, 
                   channels = 16, 
                   dropout = 0.5, 
                   l2_reg = 5e-4, 
                   learning_rate = 1e-2, 
                   es_patience = 30):
    
    train_idx,val_idx,test_idx = limit_data(labels)
    
    N = len(adata.obs.index) # number of nodes/cells
    
    # Matrix
    X = adata.X
    
    #set the mask
    train_mask = np.zeros((N,),dtype=bool)
    train_mask[train_idx] = True

    val_mask = np.zeros((N,),dtype=bool)
    val_mask[val_idx] = True

    test_mask = np.zeros((N,),dtype=bool)
    test_mask[test_idx] = True
    
    labels_encoded, classes = encode_label(labels)
    
    F = len(adata.var) # number of features (i.e., genes)
    
    # Parameters
    #channels = 16           # Number of channels in the first layer
    #dropout = 0.5           # Dropout rate for the features
    #l2_reg = 5e-4           # L2 regularization rate
    #learning_rate = 1e-2    # Learning rate
    #epochs = 200            # Number of training epochs
    #es_patience = 30        # Patience for early stopping

    # Preprocessing operations
    sc.pp.neighbors(adata, n_neighbors = neighbors, n_pcs = pcs)
    A = adata.obsp['connectivities']
    A = GCNConv.preprocess(A).astype('f4')

    # Model definition
    X_in = Input(shape=(F, ))
    fltr_in = Input((N, ), sparse=True)

    dropout_1 = Dropout(dropout)(X_in)
    graph_conv_1 = GCNConv(channels,
                             activation='relu',
                             kernel_regularizer=l2(l2_reg),
                             use_bias=False)([dropout_1, fltr_in])

    dropout_2 = Dropout(dropout)(graph_conv_1)
    graph_conv_2 = GCNConv(num_classes,
                             activation='softmax',
                             use_bias=False)([dropout_2, fltr_in])

    # Build model
    model = Model(inputs=[X_in, fltr_in], outputs=graph_conv_2)
    optimizer = Adam(lr=learning_rate)
    model.compile(optimizer=optimizer,
                  loss='categorical_crossentropy',
                  weighted_metrics=['acc'])
    model.summary()

    tbCallBack_GCN = tf.keras.callbacks.TensorBoard(
        log_dir='./Tensorboard_GCN_Treatment',
    )
    callback_GCN = [tbCallBack_GCN]

    validation_data = ([X, A], labels_encoded, val_mask)
    model.fit([X, A],
              labels_encoded,
              sample_weight=train_mask,
              epochs=epochs,
              batch_size=N,
              validation_data=validation_data,
              shuffle=False,
              callbacks=[
                  EarlyStopping(patience=es_patience,  restore_best_weights=True),
                  tbCallBack_GCN
              ])

    y_pred_all = model.predict([X, A], batch_size=N)
    
    softmax_gcn = pd.DataFrame(y_pred_all, columns = list(classes))
    labels2 = list(classes)
    
    layer_outputs = [layer.output for layer in model.layers]
    activation_model = Model(inputs=model.input, outputs=layer_outputs)
    activations = activation_model.predict([X,A],batch_size=N)

    weights = [layer.get_weights() for layer in model.layers]
    
    for group in list(classes):
        adata.obs[group] = list(softmax_gcn[group])
    return y_pred_all, softmax_gcn, activations, weights


def plot_phate(adata = None, groups = None, vmin = 0, vmax = 1):
    scv.pl.scatter(adata, basis = 'phate', color = groups, cmap = meld.utils.get_meld_cmap(), vmin = vmin, vmax = vmax)
    
def plot_weight_matrix(adata = None, weights = None, topGenes = 25, groups = 'Treatment', cmap = 'RdBu_r', dendrogram = False, figsize = (10,10), vmin = 0, vmax = 6):
    df_weights = pd.DataFrame(weights[3][0])
    df_weights.index = adata.var_names
    
    maxfeatures = []
    for gene in list(adata.var_names):
        maxfeatures.append(max(df_weights.loc[gene]))
        
    highgenes = pd.Series(maxfeatures, index=list(adata.var_names))
    
    top = highgenes.sort_values(ascending=False).iloc[0:topGenes]
    
    sc.pl.heatmap(adata, list(top.index), groupby=groups, cmap=cmap, dendrogram=dendrogram, swap_axes=True, figsize = figsize, vmin = vmin, vmax = vmax)

def Cellocluster(adata = None, activations = None, k = 3, random_state = 0):
    kmeans = KMeans(n_clusters = k, random_state = random_state).fit(activations[3])
    kmeans_encoded, classes_k = encode_label(kmeans.labels_)
    
    clusters = (pd.Series(kmeans.labels_).astype('str')).astype('category')
    
    adata.obs['Clusters'] = list(clusters)
    adata.obs['Clusters'] = adata.obs['Clusters'].astype('category')

def plotComposition(adata = None, groups = None, labels = None, stacked = True, figsize = (10,8), prediction = False):
    if (prediction == True):
        adata.obs['Prediction'] = adata.obs[list(classes)].idxmax(axis=1)
    pd.crosstab(adata.obs[groups], adata.obs[labels]).plot(kind = 'bar',
                                                          color = adata.uns[labels+'_colors'],
                                                          stacked = True, figsize = (10,8)).legend(bbox_to_anchor = (1,1.02))


def plot_violin(adata = None, labels = None, palette = 'muted'):
    labels_encoded, classes = encode_label(adata.obs[labels])
    	
    sns.set(rc={'figure.figsize':(10,8)})
    for label in classes:
    	sns.violinplot(x="Treatment", y = label, data = pd.DataFrame(adata.obs[['Treatment', 'Ctrl', 'KPT']]), order=classes,
    	palette = "muted").set(title = label + ' cells')
    	#fig = ax.get_figure()
    	
    	plt.savefig('Treatment_' + label + '_violin_plots.tif')
    	
    	plt.close()

















