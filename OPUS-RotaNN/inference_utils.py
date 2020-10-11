# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 18:41:38 2020

@author: xugang

"""

import os
import tensorflow as tf
import numpy as np
import pandas as pd

from mkinputs import PDBreader, structure, getPhiPsi, mkCSF120

ss8 = "CSTHGIEB"

ss8_dict = {}
for k,v in enumerate(ss8):
    ss8_dict[v] = k
    
def get_psp_dict():
    resname_to_psp_dict = {}
    resname_to_psp_dict['G'] = [1,4,7]
    resname_to_psp_dict['A'] = [1,3,7]
    resname_to_psp_dict['V'] = [1,7,12]
    resname_to_psp_dict['I'] = [1,3,7,12]
    resname_to_psp_dict['L'] = [1,5,7,12]
    resname_to_psp_dict['S'] = [1,2,5,7]
    resname_to_psp_dict['T'] = [1,7,15]
    resname_to_psp_dict['D'] = [1,5,7,11]
    resname_to_psp_dict['N'] = [1,5,7,14]
    resname_to_psp_dict['E'] = [1,6,7,11]
    resname_to_psp_dict['Q'] = [1,6,7,14]
    resname_to_psp_dict['K'] = [1,5,6,7,10]
    resname_to_psp_dict['R'] = [1,5,6,7,13]
    resname_to_psp_dict['C'] = [1,7,8]
    resname_to_psp_dict['M'] = [1,6,7,9]
    resname_to_psp_dict['F'] = [1,5,7,16]
    resname_to_psp_dict['Y'] = [1,2,5,7,16]
    resname_to_psp_dict['W'] = [1,5,7,18]
    resname_to_psp_dict['H'] = [1,5,7,17]
    resname_to_psp_dict['P'] = [7,19]
    return resname_to_psp_dict

def get_pc7_dict():
    resname_to_pc7_dict = {'A': [-0.350, -0.680, -0.677, -0.171, -0.170, 0.900, -0.476],
                'C': [-0.140, -0.329, -0.359, 0.508, -0.114, -0.652, 0.476],
                'D': [-0.213, -0.417, -0.281, -0.767, -0.900, -0.155, -0.635],
                'E': [-0.230, -0.241, -0.058, -0.696, -0.868, 0.900, -0.582],
                'F': [ 0.363, 0.373, 0.412, 0.646, -0.272, 0.155, 0.318],
                'G': [-0.900, -0.900, -0.900, -0.342, -0.179, -0.900, -0.900],
                'H': [ 0.384, 0.110, 0.138, -0.271, 0.195, -0.031, -0.106],
                'I': [ 0.900, -0.066, -0.009, 0.652, -0.186, 0.155, 0.688],
                'K': [-0.088, 0.066, 0.163, -0.889, 0.727, 0.279, -0.265],
                'L': [ 0.213, -0.066, -0.009, 0.596, -0.186, 0.714, -0.053],
                'M': [ 0.110, 0.066, 0.087, 0.337, -0.262, 0.652, -0.001],
                'N': [-0.213, -0.329, -0.243, -0.674, -0.075, -0.403, -0.529],
                'P': [ 0.247, -0.900, -0.294, 0.055, -0.010, -0.900, 0.106],
                'Q': [-0.230, -0.110, -0.020, -0.464, -0.276, 0.528, -0.371],
                'R': [ 0.105, 0.373, 0.466, -0.900, 0.900, 0.528, -0.371],
                'S': [-0.337, -0.637, -0.544, -0.364, -0.265, -0.466, -0.212],
                'T': [ 0.402, -0.417, -0.321, -0.199, -0.288, -0.403, 0.212],
                'V': [ 0.677, -0.285, -0.232, 0.331, -0.191, -0.031, 0.900],
                'W': [ 0.479, 0.900, 0.900, 0.900, -0.209, 0.279, 0.529],
                'Y': [ 0.363, 0.417, 0.541, 0.188, -0.274, -0.155, 0.476]}
    return resname_to_pc7_dict

resname_to_psp_dict = get_psp_dict()
resname_to_pc7_dict = get_pc7_dict()
    
def read_pssm(fname,seq):
    num_pssm_cols = 44
    pssm_col_names = [str(j) for j in range(num_pssm_cols)]
    with open(fname,'r') as f:
        tmp_pssm = pd.read_csv(f,delim_whitespace=True,names=pssm_col_names).dropna().values[:,2:22].astype(float)
    if tmp_pssm.shape[0] != len(seq):
        raise ValueError('PSSM file is in wrong format or incorrect!')
    return tmp_pssm

def read_hhm(fname,seq):
    num_hhm_cols = 22
    hhm_col_names = [str(j) for j in range(num_hhm_cols)]
    with open(fname,'r') as f:
        hhm = pd.read_csv(f,delim_whitespace=True,names=hhm_col_names)
    pos1 = (hhm['0']=='HMM').idxmax()+3
    num_cols = len(hhm.columns)
    hhm = hhm[pos1:-1].values[:,:num_hhm_cols].reshape([-1,44])
    hhm[hhm=='*']='9999'
    if hhm.shape[0] != len(seq):
        raise ValueError('HHM file is in wrong format or incorrect!')
    return hhm[:,2:-12].astype(float)

def read_ss(path):

    with open(path,'r') as f:
        ss_result = [i for i in f.readlines()]
    ss_result = ss_result[0]
    ss_len = len(ss_result)
    
    ss8 = np.zeros(ss_len*8).reshape(ss_len, 8)
    for i in range(ss_len):
        ss8[i][ss8_dict[ss_result[i]]] = 1
        
    ss3 = np.zeros((ss_len, 3))
    ss3[:,0] = np.sum(ss8[:,:3],-1)
    ss3[:,1] = np.sum(ss8[:,3:6],-1)
    ss3[:,2] = np.sum(ss8[:,6:8],-1)
    
    return ss8, ss3

def get_pp_csf120(file_path, filename, preparation_config):
    
    fasta_path = os.path.join(preparation_config["tmp_files_path"], filename +'.fasta')
    pp_path = os.path.join(preparation_config["tmp_files_path"], filename +'.pp')
    csf120_path = os.path.join(preparation_config["tmp_files_path"], filename +'.csf120')
    
    atomsData = PDBreader.readPDB(file_path) 
    residuesData = structure.getResidueData(atomsData) 
    dihedralsData = getPhiPsi.getDihedrals(residuesData)
    
    fasta = "".join([i.resname for i in residuesData])
    assert len(fasta) == len(dihedralsData)
    
    f = open(fasta_path, 'w')
    f.writelines(">" + file_path.split('/')[-1] + "\n")
    f.writelines(fasta)
    f.close()   
    
    pps = []
    for i in dihedralsData:
        pps.append([np.sin(np.deg2rad(i.pp[0])), np.cos(np.deg2rad(i.pp[0])), 
                    np.sin(np.deg2rad(i.pp[1])), np.cos(np.deg2rad(i.pp[1]))])
    f.close()
    
    np.savetxt(pp_path, np.array(pps), fmt="%.4f")
    
    mkCSF120.get_csf120(residuesData, csf120_path)

def get_pssm(filename, preparation_config):
    
    fasta_path = os.path.join(preparation_config["tmp_files_path"], filename+'.fasta')
    output_path = os.path.join(preparation_config["tmp_files_path"], filename+'.txt')
    pssm_path = os.path.join(preparation_config["tmp_files_path"], filename+'.pssm')
    
    cmd = preparation_config["psiblast_path"] + " -num_threads " + str(preparation_config["num_threads"]) + " -query " + \
            fasta_path + " -db " + preparation_config["uniref90_path"] + " -out " + \
            output_path  + " -num_iterations 3 -out_ascii_pssm " + pssm_path
            
    print (cmd)    
    
    output = os.popen(cmd).read() 
    
    if os.path.exists(output_path):
        os.remove(output_path)

def get_hhm(filename, preparation_config):
    
    fasta_path = os.path.join(preparation_config["tmp_files_path"], filename+'.fasta')
    a3m_path = os.path.join(preparation_config["tmp_files_path"], filename+'.a3m')
    hhm_path = os.path.join(preparation_config["tmp_files_path"], filename+'.hhm')
    hhr_path = os.path.join(preparation_config["tmp_files_path"], filename+'.hhr')
   
    cmd = preparation_config["hhblits_path"] + " -i " + fasta_path + \
            " -ohhm " + hhm_path + " -oa3m " + a3m_path + " -d " + preparation_config["uniclust30_path"] + \
            " -v 0 -maxres 40000 -cpu " + str(preparation_config["num_threads"]) + " -Z 0"
            
    print (cmd)    
    
    output = os.popen(cmd).read() 
    
    if os.path.exists(a3m_path):
        os.remove(a3m_path)
    if os.path.exists(hhr_path):
        os.remove(hhr_path)

def get_ss(file_path, filename, preparation_config):
    
    ss_path = os.path.join(preparation_config["tmp_files_path"], filename +'.ss')
   
    cmd = preparation_config["mkdssp_path"] + ' ' + file_path
    print (cmd) 
    
    output = os.popen(cmd).read()

    ss = []
    for i in output.split("\n"):
        if i != "" and i[0] != '#':
            ss.append(i.strip().split()[2].strip())

    f = open(ss_path, 'w')
    f.writelines("".join(ss))
    f.close()
            
def make_input(filename, preparation_config):
    """
    20pssm + 30hhm + 7pc + 19psp + 8ss + 3ss + 4pp + 120csf120 = 211
    """    
    
    fasta_path = os.path.join(preparation_config["tmp_files_path"], filename+'.fasta')
    with open(fasta_path,'r') as f:
        fasta_result = [i for i in f.readlines()]
    fasta = fasta_result[1]
    
    seq_len = len(fasta)
    
    pssm_path = os.path.join(preparation_config["tmp_files_path"], filename+'.pssm')
    hhm_path = os.path.join(preparation_config["tmp_files_path"], filename+'.hhm')
    pp_path = os.path.join(preparation_config["tmp_files_path"], filename+'.pp')
    csf120_path = os.path.join(preparation_config["tmp_files_path"], filename+'.csf120')
    ss_path = os.path.join(preparation_config["tmp_files_path"], filename+'.ss')
    input_path = os.path.join(preparation_config["tmp_files_path"], filename+'.inputs')
    
    pssm = read_pssm(pssm_path, fasta)
    hhm = read_hhm(hhm_path, fasta)
    
    ss8, ss3 = read_ss(ss_path)
    
    pc7 = np.zeros((seq_len, 7))
    for i in range(seq_len):
        pc7[i] = resname_to_pc7_dict[fasta[i]]
    
    psp = np.zeros((seq_len, 19))
    for i in range(seq_len):
        psp19 = resname_to_psp_dict[fasta[i]]
        for j in psp19:
            psp[i][j-1] = 1
    
    inputs_ = np.concatenate((pssm, hhm, pc7, psp),axis=1)
    assert inputs_.shape == (seq_len,76)
    
    pp = np.loadtxt(pp_path)
    
    csf120 = np.loadtxt(csf120_path)
    csf120 = csf120.reshape((csf120.shape[0], 10, 13))
    csf120 = csf120[:,:,1:].reshape((csf120.shape[0], 120))

    input_data = np.concatenate((inputs_, ss8, ss3, pp, csf120),-1)
    
    assert input_data.shape[-1] == 211
    
    np.savetxt(input_path, input_data, fmt="%.4f")

#=============================================================================    

def read_inputs(filenames, inputs_files_path):
    """
    20pssm + 30hhm + 7pc + 19psp + 8ss + 3ss + 4pp + 120csf120 = 211
    """
    inputs_nopadding = []
    max_len = 0
    inputs_total_len = 0
    for filename in filenames:
        inputs_ = np.loadtxt((os.path.join(inputs_files_path, filename + ".inputs")))
        
        inputs_total_len += inputs_.shape[0]
        if inputs_.shape[0] > max_len:
            max_len = inputs_.shape[0]
        inputs_nopadding.append(inputs_)
    
    inputs_padding = np.zeros(shape=(len(filenames), max_len, 211))
    inputs_mask_padding = np.ones(shape=(len(filenames), max_len))

    for i in range(len(filenames)):
        inputs_padding[i,:inputs_nopadding[i].shape[0]] = inputs_nopadding[i]
        inputs_mask_padding[i,:inputs_nopadding[i].shape[0]] = 0
        
    #(hhm - 5000) / 1000
    inputs_padding[:,:,20:50] = (inputs_padding[:,:,20:50] - 5000)/1000
        
    return inputs_padding, inputs_mask_padding, inputs_total_len

class InputReader(object):

    def __init__(self, data_list, num_batch_size, inputs_files_path):

        self.data_list = data_list
        self.inputs_files_path = inputs_files_path
        self.dataset = tf.data.Dataset.from_tensor_slices(self.data_list).batch(num_batch_size)          
        
        print ("Data Size:", len(self.data_list)) 
    
    def read_file_from_disk(self, filenames_batch):
        
        filenames_batch = [bytes.decode(i) for i in filenames_batch.numpy()]
        inputs_batch, inputs_masks_batch, inputs_total_len = \
            read_inputs(filenames_batch, self.inputs_files_path)
        
        inputs_batch = tf.convert_to_tensor(inputs_batch, dtype=tf.float32)
        inputs_masks_batch= tf.convert_to_tensor(inputs_masks_batch, dtype=tf.float32)
        
        return filenames_batch, inputs_batch, inputs_masks_batch, inputs_total_len
            
#=============================================================================    

def get_ensemble_ouput(name, predictions, x_mask, total_len):
    
    if name == "Rota":
        
        x1_predictions = []
        x2_predictions = []
        x3_predictions = []
        x4_predictions = []
        
        x1_outputs = []
        x2_outputs = []
        x3_outputs = []
        x4_outputs = []
        
        for i in predictions:
            
            # i.shape: batch, seq_len, 4
            i = i.numpy()
            
            x1_prediction = np.zeros((i.shape[0], i.shape[1], 1))
            x2_prediction = np.zeros((i.shape[0], i.shape[1], 1))
            x3_prediction = np.zeros((i.shape[0], i.shape[1], 1))
            x4_prediction = np.zeros((i.shape[0], i.shape[1], 1))

            x1_prediction[:,:,0] = np.rad2deg(np.arctan2(i[:,:,0], i[:,:,1]))
            x2_prediction[:,:,0] = np.rad2deg(np.arctan2(i[:,:,2], i[:,:,3]))
            x3_prediction[:,:,0] = np.rad2deg(np.arctan2(i[:,:,4], i[:,:,5]))
            x4_prediction[:,:,0] = np.rad2deg(np.arctan2(i[:,:,6], i[:,:,7]))
            
            x1_predictions.append(x1_prediction)
            x2_predictions.append(x2_prediction)
            x3_predictions.append(x3_prediction)
            x4_predictions.append(x4_prediction)
        
        x1_predictions = np.concatenate(x1_predictions, -1)
        x1_predictions = np.median(x1_predictions, -1)

        x2_predictions = np.concatenate(x2_predictions, -1)
        x2_predictions = np.median(x2_predictions, -1)
        
        x3_predictions = np.concatenate(x3_predictions, -1)
        x3_predictions = np.median(x3_predictions, -1)
        
        x4_predictions = np.concatenate(x4_predictions, -1)
        x4_predictions = np.median(x4_predictions, -1)
        
        x_mask = x_mask.numpy()
        max_length = x_mask.shape[1]
        for i in range(x_mask.shape[0]):
            indiv_length = int(max_length-np.sum(x_mask[i]))
            x1_outputs.append(x1_predictions[i][:indiv_length])
            x2_outputs.append(x2_predictions[i][:indiv_length])
            x3_outputs.append(x3_predictions[i][:indiv_length])
            x4_outputs.append(x4_predictions[i][:indiv_length])
        
        x1_outputs_concat = np.concatenate(x1_outputs, 0)
        x2_outputs_concat = np.concatenate(x2_outputs, 0)
        x3_outputs_concat = np.concatenate(x3_outputs, 0)
        x4_outputs_concat = np.concatenate(x4_outputs, 0)
        
        assert x1_outputs_concat.shape[0] == x2_outputs_concat.shape[0] == \
            x3_outputs_concat.shape[0] == x4_outputs_concat.shape[0] == total_len        
        
        return x1_outputs, x2_outputs, x3_outputs, x4_outputs
    
dihedral_nums ={"G":0,"A":0,"V":1,"I":2,"L":2,"S":1,"T":1,"D":2,"N":2,"E":3,"Q":3,"K":4,"R":4,"C":1,"M":3,"F":2,"Y":2,"W":2,"H":2,"P":2}
    
def output_results(filenames, x1_outputs, x2_outputs, x3_outputs, x4_outputs, preparation_config):
    
    for filename, x1_output, x2_output, x3_output, x4_output in \
        zip(filenames, x1_outputs, x2_outputs, x3_outputs, x4_outputs):
        
        fasta_path = os.path.join(preparation_config["tmp_files_path"], filename+'.fasta')
        #fasta_path = os.path.join(r'/data/xugang/opus_contact/SPOT-1D/dataset/SPOT-1D-dataset/opus_dataset/clean/fasta', 
        #                          filename+'.fasta')
        with open(fasta_path, 'r') as r:
            fasta_content = [i.strip() for i in r.readlines()]
        fasta_seq = fasta_content[1]
        
        output_path = os.path.join(preparation_config["output_path"], filename+".rota")
        f = open(output_path, 'w')
        f.write("#\tX1\tX2\tX3\tX4\n")
        
        assert x1_output.shape[0] == x2_output.shape[0] == x3_output.shape[0] == x4_output.shape[0] == len(fasta_seq)

        for idx, (res_name, x1, x2, x3, x4) in \
            enumerate(zip(fasta_seq, x1_output, x2_output, x3_output, x4_output)):
            
            null_num = 4-dihedral_nums[res_name]
            if null_num == 4:
                x1 = x2 = x3 = x4 = 182
            elif null_num == 3:
                x2 = x3 = x4 = 182
            elif null_num == 2:
                x3 = x4 = 182
            elif null_num == 1:
                x4 = 182
                
            f.write('%i\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n'%(idx+1, x1, x2, x3, x4))
            
        f.close()
    
    