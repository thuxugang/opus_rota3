# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""

from my_model import Model
import tensorflow as tf

#============================Parameters====================================
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

gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
logical_gpus = tf.config.experimental.list_logical_devices('GPU')
print(len(gpus), len(logical_gpus))
    
#============================Models====================================

model_rota1 = Model(params=params, name="rota")
model_rota1.params["save_path"] = "./models/model1"
model_rota1.load_model()  

model_rota2 = Model(params=params, name="rota")
model_rota1.params["save_path"] = "./models/model2"
model_rota2.load_model() 

model_rota3 = Model(params=params, name="rota")
model_rota1.params["save_path"] = "./models/model3"
model_rota3.load_model() 


def test_infer_step(x, x_mask):
    
    rota_predictions = []
    
    rota_prediction, _ = \
        model_rota1.inference(x, x_mask, y=None, y_mask=None, training=False)        
    rota_predictions.append(rota_prediction)

    rota_prediction, _ = \
        model_rota2.inference(x, x_mask, y=None, y_mask=None, training=False)        
    rota_predictions.append(rota_prediction)

    rota_prediction, _ = \
        model_rota3.inference(x, x_mask, y=None, y_mask=None, training=False)        
    rota_predictions.append(rota_prediction)

    return rota_predictions