# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""
import time
from my_model import Model
import tensorflow as tf
import numpy as np
from tensorflow import keras
from utils import InputReader, cal_accurarcy

if __name__ == '__main__':
    
    #parameters of training
    batch_size = 4
    epochs = 40
    early_stop = 4
    input_normalization = True
    learning_rate = 1e-3
    
    params = {}
    params["d_input"] = 211
    params["d_rota_output"] = 8
    params["dropout_rate"] = 0.4
  
    #parameters of transfomer model
    params["transfomer_layers"] = 2
    params["transfomer_num_heads"] = 1
    
    #parameters of birnn model
    params["lstm_layers"] = 4
    params["lstm_units"] = 1024
    
    params["save_path"] = r'./models'

    params["train_files_path"] = "/data/xugang/opus_rota3/github/opus_tass_datasets/Train"
    params["val_files_path"] = "/data/xugang/opus_rota3/github/opus_tass_datasets/Val"

    params["csf120_path"] = "/data/xugang/opus_rota3/github/csf120"

    gpus = tf.config.experimental.list_physical_devices('GPU')
    tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
    logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    print(len(gpus), len(logical_gpus))
    
    model_rota = Model(params=params, name="rota")

    train_list_path = "/data/xugang/opus_rota3/github/opus_tass_datasets/list_train"
    val_list_path = "/data/xugang/opus_rota3/github/opus_tass_datasets/list_val"
    
    train_reader = InputReader(data_list=train_list_path,
                               params=params,
                               num_batch_size=batch_size,
                               name="train",
                               input_norm=input_normalization, 
                               shuffle=True,
                               data_enhance=True)
    
    val_reader = InputReader(data_list=val_list_path,
                             params=params,
                             num_batch_size=batch_size,
                             name="val",
                             input_norm=input_normalization, 
                             shuffle=False,
                             data_enhance=False)
        
    lr = tf.Variable(tf.constant(learning_rate), name='lr', trainable=False)
    optimizer = keras.optimizers.Adam(lr=lr)

    def train_step(x, x_mask, y, y_mask):
        
        with tf.GradientTape() as tape:
            rota_predictions, loss = \
                model_rota.inference(x, x_mask, y, y_mask, training=True)    
            
        trainable_variables = model_rota.transformer.trainable_variables + \
            model_rota.cnn.trainable_variables + model_rota.birnn.trainable_variables
        
        gradients = tape.gradient(loss, trainable_variables)
        optimizer.apply_gradients(
            zip(gradients, trainable_variables))
        
        return loss, rota_predictions

    def infer_step(x, x_mask):
        
        rota_predictions, _ = \
            model_rota.inference(x, x_mask, y=None, y_mask=None, training=False)
            
        return rota_predictions
    
    best_val_acc = 100
    for epoch in range(epochs):
        
        #======================Train======================
        accuracy_train_x1 = []
        accuracy_train_x2 = []
        accuracy_train_x3 = []
        accuracy_train_x4 = []
        for step, filenames_batch in enumerate(train_reader.dataset):
            start_time = time.time()

            filenames, x, x_mask, y, y_mask, inputs_total_len, labels_total_len = \
                train_reader.read_file_from_disk(filenames_batch)
            
            assert inputs_total_len == labels_total_len

            loss, rota_predictions = train_step(x, x_mask, y, y_mask)
            
            mae_x1, mae_x2, mae_x3, mae_x4 = \
                cal_accurarcy("Rota", rota_predictions, y, y_mask, total_len=inputs_total_len)
                
            accuracy_train_x1.extend(mae_x1)
            accuracy_train_x2.extend(mae_x2)
            accuracy_train_x3.extend(mae_x3)
            accuracy_train_x4.extend(mae_x4)
            
            run_time = time.time() - start_time
            
            if step % 10 == 0:
                print('Epoch: %d, step: %d, loss: %3.3f, x1: %3.2f, x2: %3.2f, x3: %3.2f, x4: %3.2f, time: %3.3f'
                      % (epoch, step, loss, np.mean(accuracy_train_x1), np.mean(accuracy_train_x2), 
                          np.mean(accuracy_train_x3), np.mean(accuracy_train_x4), run_time)) 
                
        #======================Val======================
        accuracy_val_x1 = []
        accuracy_val_x2 = []
        accuracy_val_x3 = []
        accuracy_val_x4 = []
        start_time = time.time()
        for step, filenames_batch in enumerate(val_reader.dataset):
            
            filenames, x, x_mask, y, y_mask, inputs_total_len, labels_total_len = \
                val_reader.read_file_from_disk(filenames_batch)
            
            assert inputs_total_len == labels_total_len

            rota_predictions = infer_step(x, x_mask)

            mae_x1, mae_x2, mae_x3, mae_x4 = \
                cal_accurarcy("Rota", rota_predictions, y, y_mask, total_len=inputs_total_len)
                
            accuracy_val_x1.extend(mae_x1)
            accuracy_val_x2.extend(mae_x2)
            accuracy_val_x3.extend(mae_x3)
            accuracy_val_x4.extend(mae_x4)

        run_time = time.time() - start_time
        print('Epoch: %d, lr: %s, x1: %3.2f, x2: %3.2f, x3: %3.2f, x4: %3.2f, time: %3.3f'
              % (epoch, str(lr.numpy()), np.mean(accuracy_val_x1), np.mean(accuracy_val_x2), 
                  np.mean(accuracy_val_x3), np.mean(accuracy_val_x4), run_time))   
        
        if np.mean(accuracy_val_x1) < best_val_acc:
            best_val_acc = np.mean(accuracy_val_x1)
            model_rota.save_model()
        else:
            lr.assign(lr/2)
            early_stop -= 1
        
        if early_stop == 0:
            break
    
    print ("best_val_acc:", best_val_acc)
    