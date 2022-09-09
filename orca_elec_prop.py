import os
import re
import sys

import numpy as np
import pandas as pd

from cclib.parser.orcaparser import ORCA
from itertools import chain
import json

def skiplines(openfile, nlines=0):
    '''
    Function to skip nlines + 1 lines in openfile. In other words, 
    if nlines=0 it will
    go to the next line.

    Parameters
    ----------
    openfile: object.
        File object to process.
    nlines: int.
        Number of lines to skip.

    Returns
    -------
    line: string.
        Line after skipping nlines + 1 lines.
    '''

    for i in range(nlines + 1):
        line = next(openfile)

    return line

def orca_reader(infile):
    """
    Function to read an ORCA output.

    Parameters
    ----------
    infile: object.
        Orca file object to process.

    Returns
    -------
    data: cclib.Object
        cclib object parsing an Orca Output.

    Note:
    -------
    This functions needs thelibrary cclib to work.
    """
    
    obj = ORCA(log)
    data = obj.parse()

    return data

def pol_dict_to_matrix(quad):
    """
    Function 

    Parameters
    ----------
    infile: object.
  

    Returns
    -------

    """
    diag = np.eye(3, k=0)
    mat_1=np.asarray([[0,1,0],[0,0,0],[0,0,0]], dtype='float64')
    mat_2=np.asarray([[0,0,1],[0,0,0],[0,0,0]], dtype='float64')
    mat_3=np.asarray([[0,0,0],[0,0,1],[0,0,0]], dtype='float64')
  
    pol_val=[]
    for idx, values in enumerate(quad):
        a,b=values
        pol_val.append(b)
      
    pol_val=np.asarray(pol_val, dtype='float64')

    xx_yy_zz=[diag[i]*pol_val[i] for i in range(3)]
    xy_xz=mat_1*pol_val[3]+mat_2*pol_val[4]
    yz=mat_3*pol_val[5]

    return xx_yy_zz+xy_xz+yz


def quad_polarization(infile):
    '''
    Function to browse the Quadrupole moment from an Orca output file.

    Parameters
    ----------
    infile: object.
        File object to process.

    Returns
    -------
    quad: tuple arrays.
        This is a tuple with arrays of electronic, nuclear and total polarization tensor.
    '''
        
    with open(infile) as f:    
      for line in f:
         if "QUADRUPOLE MOMENT (A.U.)" in line:
             directions=skiplines(f,2)
             nuc_pol=skiplines(f,0)
             elec_pol=skiplines(f,0)
             tot_pol=skiplines(f,0)

      # Polarizability tensor directions
      dir_list=directions.split()

      elec_pos=list(zip(dir_list,elec_pol.split()[1:]))
      nuc_pos=list(zip(dir_list,nuc_pol.split()[1:]))
      tot_pos=list(zip(dir_list,tot_pol.split()[1:]))

      return pol_dict_to_matrix(elec_pos), pol_dict_to_matrix(nuc_pos), pol_dict_to_matrix(tot_pos)
      

def dipole_moment(infile):
    '''
    Function to search the dipole moment from an Orca output file.

    Parameters
    ----------
    infile: object.
        File object to process.

    Returns
    -------
    dip array: tuple.
        Returns a tuple of arrays containing electronic, nuclear and 
        total dipole moment values.
    '''
    with open(log) as f:
      for line in f:
         if "DIPOLE MOMENT" in line:
             directions=skiplines(f,1)
             dip_elec=skiplines(f,0)
             dip_nuc=skiplines(f,0)
             dip_tot=skiplines(f,1)

    dip_dir=directions.split()

    # Filling the dip dictionary
    elec_dip=dip_elec.split()[2:]
    nuc_dip=dip_nuc.split()[3:]
    tot_dip=dip_tot.split()[4:]
    
    return elec_dip, nuc_dip, tot_dip


def find_string(infile, string):
    '''
    Function to search a string in a user provided input file.

    Parameters
    ------------
    infile: object.
        File object to process. 

    string: str.
        String value to be used in the search.
    
    Returns
    -------
    dip array: list.
        Returns a list of indeces where the String has been found. 
    '''
    result=[]

    with open(log) as f:
        for idx, line in enumerate(f):
            if string in line:
              result.append(idx)

    return result

def chunk_output(infile):
    '''
    Function...

    Parameters
    ------------
    infile: object.
        File object to process. 

    string: str.
        String value to be used in the search.
    
    Returns
    -------
    dip array: list.
        Returns two lists with singlets and triplets energy levels. 
    '''

    # Searching via string in the output
    stddft_start=find_string(log,"sTD-DFT procedure...")
    stddft_end=find_string(log,"sTD-DFT done")
    stddft_triplet=find_string(log,"calculate triplets           ...   yes")

    # Creating intervales from indeces
    stddft_sections=list(zip(stddft_start,stddft_end))

    results=[np.abs(stddft_triplet[0]-stddft_start[0]),
             np.abs(stddft_triplet[0]-stddft_start[1])]
    
    results=np.array(results)
    triplet_start=np.unravel_index(np.argmin(results, axis=None), results.shape)

    return  stddft_sections[triplet_start[0]], stddft_sections[1]

def split_excitations(infile):
    '''
    Function..
    '''
    triplet_idx, singlet_idx=chunk_output(log)

    with open(infile) as fp:
      lines = fp.readlines()

    triplets=lines[triplet_idx[0]:triplet_idx[1]]
    singlets=lines[singlet_idx[0]:singlet_idx[1]]

    return triplets, singlets

def seek_string(string, array):
    """
    Function
    """
    value_start=[]
    for idx, values in enumerate(array):
       if string in values:
           value_start.append(idx)

    return value_start[0]

def split_arrays(array) :
    """
    Function to split an array and removing blank spaces.

    Parameters
    ------------
    array: list.
        User provided list. 
    
    Returns
    -------
    list: list.
        Returns a list without blank spaces and split. 
    """
    
    return [re.split("\s+", array[i].strip()) for i in range(len(array))]

def clean_arrays(array, value):
    """
    Function
   
     Parameters
    ------------
    array: list.
        User provided list.  

    value: int.
        Integer value used for slicing the list.
    
    Returns
    -------
    list: list.
        Returns a list with the fetched text between the computed indeces. 
    """
    imp_columns=len(array[value])
    valid_columns=len(array)-value

    return [array[value+i][0:imp_columns] for i in range(valid_columns-1)]

def singlet_triplet_arrays(infile):
    """
    Function to retrieve and organize the ORCA computed excitations from a valid output file. 
    
    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    list: tuple.
        Returns a tuple of lists for singlets and triplets.  
    """
    triplets, singlets = split_excitations(log)

    exc_singlets=split_arrays(singlets)
    exc_triplets=split_arrays(triplets)

    val_sing=seek_string("state", exc_singlets)
    val_trip=seek_string("state", exc_triplets)

    return clean_arrays(exc_singlets, val_sing), clean_arrays(exc_triplets, val_trip)

def array_to_df(infile):
    """
    Function to convert arrays with polarization information into Pandas 
    dataframes.
    
    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    tuple: tuple of Pandas DataFrames
        Returns a tuple of pandas dataframes holding singlets and triplets 
        results. 
    """
    singlet,triplet=singlet_triplet_arrays(infile)
    
    my_sing=singlet[1:]
    my_trip=triplet[1:]

    df_singlet = pd.DataFrame(my_sing, columns = singlet[0])
    df_triplet = pd.DataFrame(my_trip, columns = triplet[0])
    
    return df_singlet, df_triplet

def list_mo_energies(infile):
    """
    Function to list all occupied energy levels (up to HOMO).

    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    data: list
        A list storing the energy levels from a cclib object.
    """
    data=orca_reader(infile)

    return data.moenergies[0]

def last_coord(infile):
    """
    Function to retrieve the last set of coordinates from an Orca output file.
    
    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    data: list
        A list storing the last coordinates from a cclib object.
    """
    data=orca_reader(infile)
    coords=data.metadata['coords']
    natoms=data.natom
    
    return coords[natoms:]

def qm_charges(infile, method="lowdin"):
    """
    Function to retreive the lowdin charges from an Orca output file.
  
    Parameters
    ------------
    infile: object.
        File object to process. 
   
    method: string.
        String containing the method. This is lowdin or mulliken.
       

    Returns
    -------
    data: list
        A list storing the lowdin charges from a cclib object.

    """
    data=orca_reader(infile)

    # Charges from cclib
    charges=data.atomcharges
    charges_total=dict.fromkeys([str(method)], [])
    charges_total[str(method)]=np.array(charges[
        str(method)]).astype("float64").tolist()

    return charges_total

def mo_energies(infile):
    """
    Function to retrieve the occupied energy levels up to HOMO energy level.
    
    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    mo_energies: list
        A list with the energy levels for occupied orbitals from a cclib object.
    """
    data=orca_reader(infile)
  
    # Dictionary of mo_energies (WORKING)
    idx_homo=data.homos[0]
    mo_energies=dict.fromkeys(["mo_energies"], [])
    mo_energies["mo_energies"]=list_mo_energies(
        log)[0:idx_homo].tolist()

    return mo_energies

def atom3D(infile):
    """
    Function to retrieve the coordinates and elements from an Orca output file.

    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    atom3D: list
        A list with the elements and coordinates from a cclib object

    """
    data=orca_reader(infile)

    #Dictionary with atomic coordinates and elements
    atom3D=dict.fromkeys(["3D_coord"], [])
    atom3D["3D_coord"]=last_coord(log)

    return atom3D

def symbols(infile):
    """
    Function to retrieve the elements from a cclib coordinates object.

    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    symbols: dict
        A dictionary containing the elements of a cclib coordnate object.
    """
    
    data=orca_reader(log)
    coord=atom3D(log)
    
    symbols={"elements": [coord["3D_coord"][i][0]
                          for i in range(data.natom)]}

    return symbols



def energy_levels(infile):
    """
    Function to retreive and report a list of only singlet energy levels.

    Parameters
    ------------
    infile: object.
        File object to process. 
   
    Returns
    -------
    singlet_energies: list
        A list with the singlet energy levels from a pandas object.
    """
    # Function to print singlet, triplet arrays
    singlet, triplet = array_to_df(log)

    # Dictionary of singlet energies
    singlet_energies=dict.fromkeys(["singlet"], [])
    singlet_energies["singlet"]=np.array(singlet["eV"]
    ).astype("float64").tolist()

    # Dictionary of triplet energies
    triplet_energies=dict.fromkeys(["triplet"], [])
    triplet_energies["triplet"]=np.array(triplet["eV"]
    ).astype("float64").tolist()

    return singlet_energies, triplet_energies


def create_json(infile, output_name):
  """ 
   Creating dictionaries for coordinates, ground state energies
   charges, excited states (singlet and triplets), and electronic 
   properties.

   Parameters
   ------------
   infile: object.
        File object to process. 
   
   Returns
    -------
    singlet_energies: None
        A json file based on the created dictionaries. 
  """

  sin,trip=energy_levels(log)
  final = dict(chain(sin.items(), trip.items()))
  
  # Polarization matrices from inputfile 
  pol_elec, pol_nuc, pol_tot=quad_polarization(log)

  flags=["elec","nuc","tot"]
  pol=dict.fromkeys(flags,[])

  pol["elec"]=pol_elec.tolist()
  pol["nuc"]=pol_nuc.tolist()
  pol["tot"]=pol_tot.tolist()

  # Create dipole moment dictionary
  dip_elec, dip_nuc, dip_tot=dipole_moment(log)

  flags=["elec","nuc","tot"]
  dip=dict.fromkeys(flags,[])

  dip["elec"]=dip_elec
  dip["nuc"]=dip_nuc
  dip["tot"]=dip_tot
 
  
  elements=symbols(log)
  coord=atom3D(log)
  mo=mo_energies(log)
  exc_states={"excited_states": final}
  elec_prop={"dipole": dip, "polarization": pol}
  charges=qm_charges(log,"lowdin")
  
  schema={"schema_name": "orca_schema_output"}
  schema.update(elements)
  schema.update(coord)
  schema.update(mo)
  schema.update(exc_states)
  schema.update(elec_prop)
  schema.update(charges)

  with open(str(output_name)+".json","w") as f:
    json.dump(schema, f, indent=4)

if __name__ == "__main__":
    
  # Using cclib to obtain some properties
  log="output.log"
  create_json(log,"polymer")





 



 
