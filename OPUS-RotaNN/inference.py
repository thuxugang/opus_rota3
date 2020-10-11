# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""
import os
import time
from inference_utils import get_pp_csf120, get_pssm, get_hhm, get_ss, make_input, \
                            InputReader, get_ensemble_ouput, output_results
from inference_models import test_infer_step

if __name__ == '__main__':
    
    #============================Parameters====================================
    list_path = r"./list_caspfm56"
    files_path = []
    f = open(list_path)
    for i in f.readlines():
        files_path.append(i.strip())
    f.close()
   
    preparation_config = {}
    preparation_config["tmp_files_path"] = "./tmp_files"
    preparation_config["output_path"] = "./predictions"
    
    preparation_config["num_threads"] = 40
    preparation_config["psiblast_path"] = '/data/xugang/opus_contact/blast/ncbi-blast-2.10.0+/bin/psiblast'
    preparation_config["uniref90_path"] = '/data/xugang/opus_contact/uniref90/uniref90.fasta'
    preparation_config["hhblits_path"] = '/data/xugang/opus_contact/hhblits/hh-suite/build/bin/hhblits'
    preparation_config["uniclust30_path"] = '/data/xugang/opus_contact/uniclust30/uniclust30_2018_08/uniclust30_2018_08'
    
    preparation_config["mkdssp_path"] = r'./mkdssp/mkdssp'
    
    batch_size = 8
    #============================Parameters====================================
    
    
    #============================Preparation===================================
    start_time = time.time()
    filenames = []
    for file_path in files_path:
        
        filename = file_path.split('/')[-1].split('.')[0]
        
        #write pp & ci & fasta
        pp_filename = filename + '.pp'
        csf120_filename = filename + '.csf120'
        if (not os.path.exists(os.path.join(preparation_config["tmp_files_path"], pp_filename)) or
            not os.path.exists(os.path.join(preparation_config["tmp_files_path"], csf120_filename))):
            get_pp_csf120(file_path, filename, preparation_config)  
            
        pssm_filename = filename + '.pssm'
        if not os.path.exists(os.path.join(preparation_config["tmp_files_path"], pssm_filename)):
            get_pssm(filename, preparation_config)
        
        hhm_filename = filename + '.hhm'
        if not os.path.exists(os.path.join(preparation_config["tmp_files_path"], hhm_filename)):
            get_hhm(filename, preparation_config)       

        ss_filename = filename + '.ss'
        if not os.path.exists(os.path.join(preparation_config["tmp_files_path"], ss_filename)):
            get_ss(file_path, filename, preparation_config)  

        make_input(filename, preparation_config)
        
        filenames.append(filename)
        
    run_time = time.time() - start_time
    print('Preparation done..., time: %3.3f' % (run_time))  
    #============================Preparation===================================
    
    #==================================Model===================================
    start_time = time.time()
    test_reader = InputReader(data_list=filenames, 
                              num_batch_size=batch_size,
                              inputs_files_path=preparation_config["tmp_files_path"])
    
    total_lens = 0
    for step, filenames_batch in enumerate(test_reader.dataset):

        filenames, x, x_mask, inputs_total_len = \
            test_reader.read_file_from_disk(filenames_batch)
        
        total_lens += inputs_total_len
        
        rota_predictions = \
            test_infer_step(x, x_mask)
            
            
        x1_outputs, x2_outputs, x3_outputs, x4_outputs = \
            get_ensemble_ouput("Rota", rota_predictions, x_mask, inputs_total_len)    
        
        assert len(filenames) == len(x1_outputs) == len(x2_outputs) == len(x3_outputs) == len(x4_outputs)
            
        output_results(filenames, x1_outputs, x2_outputs, x3_outputs, x4_outputs, preparation_config)
        
    run_time = time.time() - start_time
    print('Prediction done..., time: %3.3f' % (run_time)) 
    #==================================Model===================================
    
    
    