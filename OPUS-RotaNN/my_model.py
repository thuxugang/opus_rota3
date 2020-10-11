# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""

import os
import tensorflow as tf
from my_transformer import Transformer, create_padding_mask
from my_rnn import BiRNN_Rota
from my_cnn import CNN
from utils import clean_inputs, compute_mse_loss

class Model(object):
    
    def __init__(self, params, name):
        
        self.params = params
        self.name = name
        
        self.transformer = Transformer(num_layers=self.params["transfomer_layers"],
                                       d_model=self.params["d_input"],
                                       num_heads=self.params["transfomer_num_heads"],
                                       rate=self.params["dropout_rate"])
        
        self.cnn = CNN()  
        
        self.birnn = BiRNN_Rota(num_layers=self.params["lstm_layers"],
                              units=self.params["lstm_units"],
                              rate=self.params["dropout_rate"],
                              rota_output=self.params["d_rota_output"])
        print ("use rota model...")

    def inference(self, x, x_mask, y, y_mask, training):

        encoder_padding_mask = create_padding_mask(x_mask)    
        
        x = clean_inputs(x, x_mask, self.params["d_input"])
        
        transformer_out = self.transformer(x, encoder_padding_mask, training=training)
        cnn_out = self.cnn(x, training=training)
        x = tf.concat((x, cnn_out, transformer_out), -1)
        
        x = clean_inputs(x, x_mask, 3*self.params["d_input"])
        
        rota_predictions = \
            self.birnn(x, x_mask, training=training) 
        loss = None
        if training == True:
            loss = compute_mse_loss(rota_predictions, y[:,:,15:23], y_mask[:,:,15:23])
        
        return rota_predictions, loss

    def save_model(self):
        print ("save model:", self.name)
        self.transformer.save_weights(os.path.join(self.params["save_path"], self.name + '_trans_model_weight'))
        self.cnn.save_weights(os.path.join(self.params["save_path"], self.name + '_cnn_model_weight'))
        self.birnn.save_weights(os.path.join(self.params["save_path"], self.name + '_birnn_model_weight'))

    def load_model(self):
        print ("load model:", self.name)
        self.transformer.load_weights(os.path.join(self.params["save_path"], self.name + '_trans_model_weight'))
        self.cnn.load_weights(os.path.join(self.params["save_path"], self.name + '_cnn_model_weight'))
        self.birnn.load_weights(os.path.join(self.params["save_path"], self.name + '_birnn_model_weight'))





