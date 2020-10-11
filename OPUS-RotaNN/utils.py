# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""

import os
import tensorflow as tf
import numpy as np
from tensorflow import keras
import scipy.stats as stats

def read_filenames(data_list):
    
    filenames = []
    f = open(data_list, 'r')
    for i in f.readlines():
        if i.strip() != "":
            filenames.append(i.strip())
    f.close()

    return filenames

ratio = 0.25
def get_enhancement(inputs, index=None):
    
    if index != None:
        return inputs[index[0]:index[1]]
    else:
        length = inputs.shape[0]
        # about half
        if np.random.randint(0,2) == 0:
            return inputs, [0, length]
        else:
            start = np.random.randint(0, int(length*ratio))
            end = length - np.random.randint(0, int(length*ratio))       
            return inputs[start:end], [start, end]

def read_inputs(filenames, params, name, data_enhance, input_norm):
    """
    .inputs
    20pssm + 30hhm + 7pc + 19psp
    .labels
    8ss(one-hot) + 3csf(double) + [2*(phi+psi) + 2*(x1+x2+x3+x4)](sin,cos) + asa + real_phipsidihedrals
    8 + 3 + 4 + 8 + 1 + 6
    
    inputs_
    20pssm + 30hhm + 7pc + 19psp + ss8 + ss3 + phipsi4 + csf120
    """
    inputs_nopadding = []
    max_len = 0
    inputs_total_len = 0
    indices = []
    for filename in filenames:
        
        if name == 'train':
            dir_name = params["train_files_path"]
        elif name == 'val':
            dir_name = params["val_files_path"]
               
        inputs_ = np.loadtxt(os.path.join(dir_name, filename + ".inputs"))
        
        labels_ = np.loadtxt(os.path.join(dir_name, filename + ".labels"))
        labels_ss3 = np.zeros((labels_.shape[0], 3))
        labels_ss3[:,0] = np.sum(labels_[:,:3],-1)
        labels_ss3[:,1] = np.sum(labels_[:,3:6],-1)
        labels_ss3[:,2] = np.sum(labels_[:,6:8],-1)
        labels_ = np.concatenate((labels_, labels_ss3),-1)
        
        csf120_ = np.loadtxt((os.path.join(params["csf120_path"], filename + ".csf120")))
        csf120_ = csf120_.reshape((inputs_.shape[0], 10, 13))
        csf120_ = csf120_[:,:,1:].reshape((inputs_.shape[0], 120))

        inputs_ = np.concatenate((inputs_, labels_[:,:8], labels_[:,30:33], labels_[:,11:15], csf120_),-1)
        assert inputs_.shape == (inputs_.shape[0], params["d_input"])
        
        if data_enhance:
            inputs_, index = get_enhancement(inputs_)
            indices.append(index)
            
        inputs_total_len += inputs_.shape[0]
        if inputs_.shape[0] > max_len:
            max_len = inputs_.shape[0]
        inputs_nopadding.append(inputs_)
    
    inputs_padding = np.zeros(shape=(len(filenames), max_len, params["d_input"]))
    inputs_mask_padding = np.ones(shape=(len(filenames), max_len))

    for i in range(len(filenames)):
        inputs_padding[i,:inputs_nopadding[i].shape[0]] = inputs_nopadding[i]
        inputs_mask_padding[i,:inputs_nopadding[i].shape[0]] = 0
        
    if input_norm:
        #(hhm - 5000) / 1000
        inputs_padding[:,:,20:50] = (inputs_padding[:,:,20:50] - 5000)/1000
        
    return inputs_padding, inputs_mask_padding, inputs_total_len, indices

def read_labels(filenames, params, name, data_enhance, indices):
    """
    ss_labels = labels[:,:,:8]
    csf_labels = labels[:,:,8:11]
    phipsi_labels = labels[:,:,11:15]
    dihedrals_labels = labels[:,:,15:23]
    asa_labels = labels[:,:,23]
    real_phipsidihedrals=labels[:,24:30]
    """
    labels_nopadding = []
    masks_nopadding = []
    max_len = 0
    labels_total_len = 0
    for idx, filename in enumerate(filenames):

        if name == 'train':
            dir_name = params["train_files_path"]
        elif name == 'val':
            dir_name = params["val_files_path"]
        
        labels_ = np.loadtxt(os.path.join(dir_name, filename + ".labels"))
        masks_ = np.loadtxt(os.path.join(dir_name, filename + ".labels_mask"))
        
        if data_enhance:
            labels_ = get_enhancement(labels_, indices[idx])
            masks_ = get_enhancement(masks_, indices[idx])
            
        assert labels_.shape[0] == masks_.shape[0]
        labels_total_len += labels_.shape[0]
        if labels_.shape[0] > max_len:
            max_len = labels_.shape[0]
        labels_nopadding.append(labels_)
        masks_nopadding.append(masks_)
        
    labels_padding = np.zeros(shape=(len(filenames), max_len, 30))
    masks_padding = np.ones(shape=(len(filenames), max_len, 30))

    for i in range(len(filenames)):
        labels_padding[i,:labels_nopadding[i].shape[0]] = labels_nopadding[i]
        masks_padding[i,:masks_nopadding[i].shape[0]] = masks_nopadding[i]
        
    return labels_padding, masks_padding, labels_total_len

class InputReader(object):

    def __init__(self, data_list, params, \
                 num_batch_size, name, input_norm=False, shuffle=False, data_enhance=False):

        self.filenames = read_filenames(data_list)
        self.params = params
        self.name = name
        self.input_norm = input_norm
        self.data_enhance = data_enhance
        
        if self.data_enhance:
            print ("use data enhancement...")
            
        if shuffle:
            self.dataset = tf.data.Dataset.from_tensor_slices(self.filenames) \
                .shuffle(len(self.filenames)).batch(num_batch_size)
        else:
             self.dataset = tf.data.Dataset.from_tensor_slices(self.filenames) \
                .batch(num_batch_size)          
        
        print ("Data Size:", len(self.filenames)) 
    
    def read_file_from_disk(self, filenames_batch):
        
        filenames_batch = [bytes.decode(i) for i in filenames_batch.numpy()]
        inputs_batch, inputs_masks_batch, inputs_total_len, indices = \
            read_inputs(filenames_batch, self.params, self.name, self.data_enhance, self.input_norm)
        labels_batch, labels_masks_batch, labels_total_len = \
            read_labels(filenames_batch, self.params, self.name, self.data_enhance, indices) 
        
        inputs_batch = tf.convert_to_tensor(inputs_batch, dtype=tf.float32)
        inputs_masks_batch= tf.convert_to_tensor(inputs_masks_batch, dtype=tf.float32)
        labels_batch = tf.convert_to_tensor(labels_batch, dtype=tf.float32)
        labels_masks_batch= tf.convert_to_tensor(labels_masks_batch, dtype=tf.float32)
        
        return filenames_batch, inputs_batch, inputs_masks_batch, \
            labels_batch, labels_masks_batch, inputs_total_len, labels_total_len
    
mse_loss_func = keras.losses.MeanSquaredError()

def compute_mse_loss(predictions, labels, labels_mask):
    
    # labels.shape: batch, seq_len, 4
    # labels_mask.shape: batch, seq_len, 4
    # predictions.shape: batch, seq_len, 4

    labels = tf.reshape(labels, (tf.shape(labels)[0]*tf.shape(labels)[1]*tf.shape(labels)[2],))
    labels_mask = tf.reshape(labels_mask, (tf.shape(labels_mask)[0]*tf.shape(labels_mask)[1]*tf.shape(labels_mask)[2],))
    predictions = tf.reshape(predictions, (tf.shape(predictions)[0]*tf.shape(predictions)[1]*tf.shape(predictions)[2],))

    # labels_mask.shape: batch, seq_len
    indices = tf.squeeze(tf.where(tf.math.equal(labels_mask, 0)), 1)
    labels_ = tf.gather(labels, indices)
    predictions_ = tf.gather(predictions, indices)

    # loss_.shape: batch*seq_len*3
    loss_ = mse_loss_func(labels_, predictions_)
    
    return loss_  

def cal_accurarcy(name, predictions, labels, labels_mask, total_len):

    if name == "Rota":

        predictions = tf.reshape(predictions, (tf.shape(predictions)[0]*tf.shape(predictions)[1], 8))
        
        x1_labels = labels[:,:,26]
        x1_labels_mask = labels_mask[:,:,26]
        x2_labels = labels[:,:,27]
        x2_labels_mask = labels_mask[:,:,27]
        x3_labels = labels[:,:,28]
        x3_labels_mask = labels_mask[:,:,28]
        x4_labels = labels[:,:,29]
        x4_labels_mask = labels_mask[:,:,29]
        
        x1_labels = tf.reshape(x1_labels, (tf.shape(labels)[0]*tf.shape(labels)[1],))
        x1_labels_mask = tf.reshape(x1_labels_mask, (tf.shape(labels_mask)[0]*tf.shape(labels_mask)[1],))
        x2_labels = tf.reshape(x2_labels, (tf.shape(labels)[0]*tf.shape(labels)[1],))
        x2_labels_mask = tf.reshape(x2_labels_mask, (tf.shape(labels_mask)[0]*tf.shape(labels_mask)[1],))
        x3_labels = tf.reshape(x3_labels, (tf.shape(labels)[0]*tf.shape(labels)[1],))
        x3_labels_mask = tf.reshape(x3_labels_mask, (tf.shape(labels_mask)[0]*tf.shape(labels_mask)[1],))
        x4_labels = tf.reshape(x4_labels, (tf.shape(labels)[0]*tf.shape(labels)[1],))
        x4_labels_mask = tf.reshape(x4_labels_mask, (tf.shape(labels_mask)[0]*tf.shape(labels_mask)[1],))        
        
        indices = tf.squeeze(tf.where(tf.math.equal(x1_labels_mask, 0)), 1)
        x1_labels_ = tf.gather(x1_labels, indices)
        predictions_ = tf.gather(predictions, indices)
        assert x1_labels_.shape[0] == predictions_.shape[0]
        
        x1_labels_ = x1_labels_.numpy()
        predictions_ = predictions_.numpy()

        # predictions.shape: batch*seq_len, 2
        predictions_ = np.rad2deg(np.arctan2(predictions_[:,0], predictions_[:,1]))
        x1_diff = x1_labels_ - predictions_
        x1_diff[np.where(x1_diff<-180)] += 360
        x1_diff[np.where(x1_diff>180)] -= 360
        mae_x1 = np.abs(x1_diff)
        
        indices = tf.squeeze(tf.where(tf.math.equal(x2_labels_mask, 0)), 1)
        x2_labels_ = tf.gather(x2_labels, indices)
        predictions_ = tf.gather(predictions, indices)
        assert x2_labels_.shape[0] == predictions_.shape[0]
        
        x2_labels_ = x2_labels_.numpy()
        predictions_ = predictions_.numpy()

        # predictions.shape: batch*seq_len, 2
        predictions_ = np.rad2deg(np.arctan2(predictions_[:,2], predictions_[:,3]))
        
        x2_diff = x2_labels_ - predictions_
        x2_diff[np.where(x2_diff<-180)] += 360
        x2_diff[np.where(x2_diff>180)] -= 360
        mae_x2 = np.abs(x2_diff)

        indices = tf.squeeze(tf.where(tf.math.equal(x3_labels_mask, 0)), 1)
        x3_labels_ = tf.gather(x3_labels, indices)
        predictions_ = tf.gather(predictions, indices)
        assert x3_labels_.shape[0] == predictions_.shape[0]
        
        x3_labels_ = x3_labels_.numpy()
        predictions_ = predictions_.numpy()

        # predictions.shape: batch*seq_len, 2
        predictions_ = np.rad2deg(np.arctan2(predictions_[:,4], predictions_[:,5]))
        
        x3_diff = x3_labels_ - predictions_
        x3_diff[np.where(x3_diff<-180)] += 360
        x3_diff[np.where(x3_diff>180)] -= 360
        mae_x3 = np.abs(x3_diff)

        indices = tf.squeeze(tf.where(tf.math.equal(x4_labels_mask, 0)), 1)
        x4_labels_ = tf.gather(x4_labels, indices)
        predictions_ = tf.gather(predictions, indices)
        assert x4_labels_.shape[0] == predictions_.shape[0]
        
        x4_labels_ = x4_labels_.numpy()
        predictions_ = predictions_.numpy()

        # predictions.shape: batch*seq_len, 2
        predictions_ = np.rad2deg(np.arctan2(predictions_[:,6], predictions_[:,7]))
        
        x4_diff = x4_labels_ - predictions_
        x4_diff[np.where(x4_diff<-180)] += 360
        x4_diff[np.where(x4_diff>180)] -= 360
        mae_x4 = np.abs(x4_diff)
        
        return mae_x1, mae_x2, mae_x3, mae_x4

def clean_inputs(x, x_mask, dim_input):
    # set 0
    # x.shape: batch, seq_len, dim_input
    # x_mask.shape: batch, seq_len
    x_mask = tf.tile(x_mask[:,:,tf.newaxis], [1, 1, dim_input])
    x_clean = tf.where(tf.math.equal(x_mask, 0), x, x_mask-1)
    return x_clean
