# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""

import tensorflow as tf
from tensorflow import keras


class CNN(keras.Model):
    def __init__(self):
        
        super(CNN, self).__init__()
        
        self.cnn5x5_layer1 = keras.layers.Conv2D(16, kernel_size=(5,5), padding='SAME')
        self.cnn5x5_layer2 = keras.layers.Conv2D(32, kernel_size=(5,5), padding='SAME')
        self.cnn5x5_layer3 = keras.layers.Conv2D(64, kernel_size=(5,5), padding='SAME')

        self.cnn3x3_layer1 = keras.layers.Conv2D(32, kernel_size=(3,3), padding='SAME')
        self.cnn3x3_layer2 = keras.layers.Conv2D(16, kernel_size=(3,3), padding='SAME')
        self.cnn3x3_layer3 = keras.layers.Conv2D(8, kernel_size=(3,3), padding='SAME')
        self.cnn3x3_layer4 = keras.layers.Conv2D(1, kernel_size=(3,3), padding='SAME')
        
        self.bn_layer1 = keras.layers.BatchNormalization()
        self.bn_layer2 = keras.layers.BatchNormalization()
        self.bn_layer3 = keras.layers.BatchNormalization()
        
    def call(self, x, training):
        
        x = tf.expand_dims(x, -1)
        
        #batch, seq, 94, 11
        x = self.cnn5x5_layer1(x)
        x = self.cnn5x5_layer2(x)
        x = self.bn_layer1(x, training=training)
        x = tf.nn.relu(x)

        x = self.cnn5x5_layer3(x)
        x = self.cnn3x3_layer1(x)
        x = self.bn_layer2(x, training=training)
        x = tf.nn.relu(x)

        x = self.cnn3x3_layer2(x)
        x = self.cnn3x3_layer3(x)
        x = self.bn_layer3(x, training=training)
        x = tf.nn.relu(x)
        
        #batch, seq, 94, 1
        x = self.cnn3x3_layer4(x)
        cnn_output = tf.reduce_mean(x, axis=-1)
        
        return cnn_output
