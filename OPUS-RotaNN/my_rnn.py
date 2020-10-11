# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import LSTM, Bidirectional

class BiRNN_Rota(keras.Model):
    def __init__(self, num_layers, units, rate, 
                 rota_output):
        super(BiRNN_Rota, self).__init__()
        
        self.num_layers = num_layers
        
        self.birnn_layers = [
            Bidirectional(LSTM(units, dropout=rate, return_sequences=True), merge_mode='concat')
            for _ in range(self.num_layers)]
        
        self.dropout = keras.layers.Dropout(rate)
        
        self.rota_layer = keras.layers.Dense(rota_output)
        
    def call(self, x, x_mask, training):
        
        # x.shape (batch_size, max_seq_length, embeded_size)
        # x_mask.shape (batch_size, timesteps).
        x_mask = tf.math.logical_not(tf.math.equal(x_mask, 1))
        
        for i in range(self.num_layers):
            x = self.birnn_layers[i](x, mask=x_mask, training=training)
        
        x = tf.nn.relu(x)
        x = self.dropout(x, training=training)
        
        # predictions.shape: (batch_size, input_seq_len, output_size)
        rota_predictions = self.rota_layer(x)

        return rota_predictions
