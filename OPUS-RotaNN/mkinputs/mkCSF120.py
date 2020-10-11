# -*- coding: utf-8 -*-
"""
Created on Wed May 11 08:43:50 2016

@author: Xu Gang
"""

from __future__ import division

import os
import numpy as np
from mkinputs import PDBreader

def getConservedInfo(residuesData):
    
    csf_data = []
    res_len = len(residuesData)
    for idx in range(res_len):
        ref_id = idx + 3
        if ref_id >= res_len:
            csf_data.append(np.array([0,0,0]))
        else:
            csf_data.append(transCoordinate(residuesData[ref_id].atoms["CA"].position, \
                residuesData[ref_id].atoms["C"].position, residuesData[ref_id].atoms["O"].position, \
                residuesData[idx].atoms["C"].position))
       
    return np.around(csf_data, decimals=2)
        
def transCoordinate(atom_ca_ref, atom_c_ref, atom_o_ref, atom_c):
    
    ref = atom_ca_ref
    c_ref_new = atom_c_ref - ref
    o_ref_new = atom_o_ref - ref
    c_new = atom_c - ref
  
    #c-ca
    x_axis = c_ref_new/np.linalg.norm(c_ref_new)
    
    c_o = o_ref_new - c_ref_new
    
    y_axis = c_o - (x_axis.dot(c_o)/x_axis.dot(x_axis) * x_axis)
    y_axis = y_axis/np.linalg.norm(y_axis)

    z_axis = np.cross(x_axis,y_axis)
    
    rotation_matrix = np.array([x_axis[0],y_axis[0],z_axis[0],0,x_axis[1],y_axis[1],z_axis[1],0,x_axis[2],y_axis[2],z_axis[2],0,0,0,0,1]).reshape(4,4)

    new = np.array([c_new[0],c_new[1],c_new[2],1]).dot(rotation_matrix)

    return np.array([new[0],new[1],new[2]])

def takeSecond(elem):
    return elem[1]

def get_contactlist(residuesData):
    
    contactlist = []
    length = len(residuesData)
    for i in range(length):
        residue_a = residuesData[i]
        a_ca = residue_a.atoms["CA"].position
        
        dis_lists = []
        for j in range(length):
            if i == j:
                continue
            residue_b = residuesData[j]
            b_ca = residue_b.atoms["CA"].position
            
            aa_distance = np.linalg.norm(a_ca - b_ca)
            dis_lists.append([j, aa_distance])
        
        assert len(dis_lists) == length - 1
        dis_lists.sort(key=takeSecond)
        
        dis_lists = dis_lists[:10]
        dis_csf_lists = []
        for item in dis_lists:

            n_data = transCoordinate(residuesData[i].atoms["CA"].position, \
                residuesData[i].atoms["C"].position, residuesData[i].atoms["O"].position, \
                residuesData[item[0]].atoms["N"].position)
            n_data = np.around(n_data, decimals=2)  
            
            ca_data = transCoordinate(residuesData[i].atoms["CA"].position, \
                residuesData[i].atoms["C"].position, residuesData[i].atoms["O"].position, \
                residuesData[item[0]].atoms["CA"].position)
            ca_data = np.around(ca_data, decimals=2)  

            c_data = transCoordinate(residuesData[i].atoms["CA"].position, \
                residuesData[i].atoms["C"].position, residuesData[i].atoms["O"].position, \
                residuesData[item[0]].atoms["C"].position)
            c_data = np.around(c_data, decimals=2)  

            o_data = transCoordinate(residuesData[i].atoms["CA"].position, \
                residuesData[i].atoms["C"].position, residuesData[i].atoms["O"].position, \
                residuesData[item[0]].atoms["O"].position)
            o_data = np.around(o_data, decimals=2)  
            
            dis_csf_lists.extend(np.concatenate(([item[0]], n_data, ca_data, c_data, o_data)))
        
        contactlist.append(dis_csf_lists)
        
    contactlist = np.array(contactlist)
    assert contactlist.shape == (length, 130)
    
    return contactlist
    
def get_csf120(residuesData, csf120_path):
        
    contactlist = get_contactlist(residuesData)
    np.savetxt(csf120_path, contactlist, fmt="%.2f")
        


    