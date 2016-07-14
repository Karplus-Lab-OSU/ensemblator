#!/usr/bin/env python


from Tkinter import *
import tkFileDialog
import tkMessageBox
import io
import ntpath
import Bio
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.SVDSuperimposer import SVDSuperimposer
import os
import sys
import numpy as np
import time
import itertools
from itertools import combinations
from itertools import permutations
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from operator import itemgetter
import platform
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
import random
from scipy.cluster.vq import kmeans, vq
from sklearn.cluster import AgglomerativeClustering

os.path.realpath(os.path.dirname(sys.argv[0]))



class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state="normal")
        self.widget.insert("end", str, (self.tag,))
        self.widget.see(END)
        self.widget.configure(state="disabled")
        self.widget.update_idletasks()

# definition for a type of dictionary that will allow new dictionaries to
# be added on demand
class NestedDict(dict):
    def __missing__(self, key):
        self[key] = NestedDict()
        return self[key]


# function for iterating
def grouped(iterable, n):
    # s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1),
    # (s2n,s2n+1,s2n+2,...s3n-1), ...
    return izip(*[iter(iterable)]*n)



class options:
    def __init__(self):
        self.input = []
        self.outout = ""
        self.pwd = ""
        self.permissive = False
        self.semipermissive = 0
        self.align = False
        self.template = ""
        self.chain = "X"
        self.model = "0"
        self.percent = 0.7
        self.dcut = 2.5
        self.auto = False
        self.maxclust = 6
        self.groupm = StrVar()
        self.groupn = StrVar()
        self.avg = False
        self.color = False




################################################################################
#### Functions ####
################################################################################



def aligner(pdb):
    pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)
    structure = pdb_reader.get_structure("temp", pdb)

    ref_model = structure[0]
    for alt_model in structure :
        try:
            #Get list of all atoms:
            ref_atoms = []
            alt_atoms = []
            for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
                for ref_res, alt_res in zip(ref_chain, alt_chain) :
                    for atom in ref_res:
                        ref_atoms.append(ref_res[atom.id])
                    for atom in alt_res:
                        alt_atoms.append(alt_res[atom.id])

            #Align these paired atom lists:
            super_imposer = Superimposer()
            super_imposer.set_atoms(ref_atoms, alt_atoms)

            if ref_model.id == alt_model.id :
                #Check for self/self get zero RMS, zero translation
                #and identity matrix for the rotation.
                assert \
                np.abs(super_imposer.rms) < 0.0000001
                assert \
                np.max(np.abs(super_imposer.rotran[1])) < 0.000001
                assert \
                np.max(np.abs(super_imposer.rotran[0]) - \
                                            np.identity(3)) < 0.000001
            else :
                #Update the structure by moving all the atoms in
                #this model (not just the ones used for the alignment)
                super_imposer.apply(alt_model.get_atoms())

            print "RMS(first model, model %i) = %0.2f" \
                    % (alt_model.id, super_imposer.rms)
        except:
            print("Failed to align model " + str(alt_model.id) + "." +
                  " Consider removal to improve final analysis."
                  )

    io=Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb)

    # fix the pdb file format
    filein = open(pdb, 'r')
    temp = []
    for line in filein:
        if line[0:3] == 'END' and line[0:6] != 'ENDMDL':
            pass
        else:
            temp.append(line)
    temp.append("END")
    filein.close()
    os.remove(pdb)
    fileout = open(pdb, 'w')
    for line in temp:
        fileout.write(line)
    fileout.close()


# eeprep function, is passed a list of xray-prepped files
# It's purpose is to check for equvilence of each atom type at each
# residue in all the structures
# it removes non equivelent atoms from the ensemble
def eeprep(pdbs, bad_files):

    pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)
    residues = {}
    remove_residues = {}
    atoms = {}
    structures = {}
    counter = 0
    legend_dict = {}


    for pdb in pdbs:
        # read in structure using biopython
        structure = pdb_reader.get_structure("temp", pdb)

        #create a dictionary that stores each structure seperatly with a
        # unique ID
        structures[counter] = structure
        legend_dict[counter] = pdb
        counter = counter + 1

        for residue in structure[0]['A']:
            # lists for atoms, and atom types
            atom_list_full = list()
            atom_name_list = list()
            atom_remove_list = list()
            # get just the residue number
            resnum = residue.get_id()[1]
            # if that residue is aleady in the residue dictionary,
            # append a list of all the atom types in this residue in
            # this chain, in this model, in this structure
            if resnum in residues:
                atom_list_full.append(residue.get_unpacked_list())
                for atom in atom_list_full[0]:
                    if atom.element == 'H':
                        atom_remove_list.append(atom.id)
                    else:
                        atom_name_list.append(atom.id)
                residues[resnum].append(set(atom_name_list))
                remove_residues[resnum].append(atom_remove_list)
            else:
                atom_list_full.append(residue.get_unpacked_list())
                for atom in atom_list_full[0]:
                    if atom.element == 'H':
                        atom_remove_list.append(atom.id)
                    else:
                        atom_name_list.append(atom.id)
                residues[resnum] = list()
                residues[resnum].append(set(atom_name_list))
                remove_residues[resnum] = list()
                remove_residues[resnum].append(atom_remove_list)

    # permutes through each set of atom types and identifies not common atoms
    # by subtracting sets to get sets of unique atoms

    # check all structures have the same residues
    resnum_count_list = list()
    for resnum in residues:
        resnum_count_list.append(len(residues[resnum]))
    try:
        n = max(resnum_count_list)
    except:
        for filename in bad_files:
            os.remove(filename)
        print("\n\n\nFATAL ERROR: \n"
              "It is likely that there are no files without gaps."
              " Please rerun while allowing some or all models with gaps"
              " . Especially consider doing an alignment. \n\n\n"
              )

    new_residues = dict()

    # append all atoms to a list if less than all members of ensemble have
    # that residue
    for resnum in residues:
        if len(residues[resnum]) < n:
            new_residues[resnum] = residues[resnum]

    removal_dict = {}
    for resnum in residues:
        # removes atoms when all structures have that residue
        for x,y in permutations(residues[resnum], r = 2):
            if resnum in removal_dict:
                # extend here, to just have one list
                # this set(x-y) is a set of any atoms in one structure
                # but not the other
                removal_dict[resnum].extend(list(set(x-y)))
            else:
                removal_dict[resnum] = list()
                removal_dict[resnum].extend(list(set(x-y)))

    # removes all atoms at this residue when some structures are missing
    # this residue
    for resnum in new_residues:
        if resnum in removal_dict:
            for atom in new_residues[resnum]:
                # extend here, to just have one list
                removal_dict[resnum].extend(list(atom))
        else:
            for atom in new_residues[resnum]:
                # extend here, to just have one list
                removal_dict[resnum] = list(atom)
    # removes all hydrogens
    for resnum in remove_residues:
        if resnum in removal_dict:
            for atom in remove_residues[resnum]:
                # extend here, to just have one list
                removal_dict[resnum].extend(list(atom))
        else:
            for atom in remove_residues[resnum]:
                removal_dict[resnum] = list(atom)

    # removes duplicate unique atoms
    for resnum in removal_dict:
        removal_dict[resnum] = set(removal_dict[resnum])

    # actual removal occurs here
    for resnum in removal_dict:
        for atom in removal_dict[resnum]:
            for key in structures:
                structure = structures[key]
                resnum = int(resnum)
                atom = atom

                #remove this atom
                try:
                    # remove this atom
                    structure[0]["A"][resnum].detach_child(atom)
                except:
                    pass


    # start building the structure
    io = PDBIO()


    # get a list of keys from the structures dictionary, and sort it
    # alphabetically, so that the final product will be more sensible

    sorted_structures = []
    for key in structures:
        # get a list of the filenames, without the "prepped_" in front
        sorted_structures.append(legend_dict[key][8:len(legend_dict[key])])
    # get the index order to call these in alphabetical order, put those
    # in a list (so now I can call the dictionaries in this order, using
    # these index numbers as the keys.)
    sorted_keys = sorted(range(len(sorted_structures)), \
                                            key=lambda k: sorted_structures[k])

    # need a counter that is just based on order, to use here, so that model 0
    # will be the first occuring one, etc.
    order_counter = 0
    # this is needed so that we can access the actual index while iterating
    # over the sorted list
    # also need the list in order to iterate over a non-random order of the
    # indicies
    # I'm aware that this is embarassingly cumbersome
    order_dict = {}
    order_list = []
    for key in sorted_keys:
        io.set_structure(structures[key])
        io.save(str(order_counter)+'_out.pdb')
        order_dict[order_counter] = key
        order_list.append(order_counter)
        order_counter += 1

    outputname = options.output

    outfile = open(outputname, 'w')
    outfile.write("")
    outfile.close

    # formatting the output files here
    counter = 0
    for key in order_dict:
        filename = str(key)+'_out.pdb'
        if backbone_scan(filename) == True:
            os.remove(filename)
            # indicate which conformations were bad
            legend_dict[order_dict[key]] = 'REMOVED FROM FINAL ENSEMBLE'
        else:
            infile = open(filename,'r')
            temp = []
            # formatting to ensure that the MODEL line is correct
            modeltag = str(counter)
            while len(modeltag) < 9:
                modeltag = " " + modeltag
            if len(modeltag) == 9:
                modeltag = "MODEL" + modeltag + "\n"
            # write model line
            temp.append(modeltag)
            # write all the atom lines
            # sort them so that they will never be out of order compared to each other
            all_atom_lines = {}
            all_atom_deets = []
            for line in infile:
                if line[0:6] == 'ATOM  ':
                    atom_deets = str(str(line[22:26].replace(" ", "")) + \
                                     str(line[12:16].replace(" ", "")))
                    all_atom_deets.append(atom_deets)
                    all_atom_lines[atom_deets] = line
            infile.close()
            all_atom_deets.sort()
            for deets in all_atom_deets:
                temp.append(all_atom_lines[deets])
            # write endmdl line
            temp.append("ENDMDL\n")
            os.remove(filename)
            outfile = open(outputname, 'a')
            for line in temp:
                outfile.write(line)
            outfile.close()
            counter += 1

    # need to cap the file with 'END'
    outfile = open(outputname, 'a')
    outfile.write("END   \n")
    outfile.close()

    # rewrites the files to correctly format them, removing the alt_conf id
    temp = []
    with open(outputname, 'r') as infile:

        for line in infile:
            if line[0:6] == 'ATOM  ':
                line = line[0:16] + " " + line[17:len(line)]
                temp.append(line)
            else:
                temp.append(line)

    os.remove(outputname)
    with open(outputname, 'w') as outfile:
        for line in temp:
            outfile.write(line)


    # now sort the legend_dict to reflect the new order of _out files
    sorted_legend_dict = {}
    for key in order_list:
        sorted_legend_dict[key] = legend_dict[order_dict[key]]

    return sorted_legend_dict

# xray-prep function, reads in a single pdb file at a time,
# also needs an output name to iterate over
# many variable names in Esperanto, sorry
def xrayprep(pdb, output):

    pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)

    # reads the pdb file
    try:
        strukturo = pdb_reader.get_structure("temp", pdb)
    except:
        print("\n\n\nFailed to read file:"
                 " '" + str(pdb) + "'. "
                 "This means there is something wrong with it."
                 " Often, using only the ATOM lines from this file will"
                 " fix this problem. If there are multiple models, be"
                 " careful getting just the ATOM lines. You will"
                 " want to get the MODEL lines as well.\n\n\n"
                 "\n\n\nIMPORTANT: The most common cause for this problem"
                 " is that the 'non-amino acid' atoms, ie. 'Cl','Ca','Se',"
                 " are formatted as seen previously. To be correct they"
                 " need to be written in all caps: 'CL','CA','SE'.\n\n\n"
                )

    cxenaro = []
    alitiparo = ['A']
    modelaro = []
    # get the list of all alternate conformation types
    for modelo in strukturo:
        modelaro.append(modelo.id)
        for cxeno in modelo:
            cxenaro.append(cxeno)
            for aminacido in cxeno:
                for atomo in aminacido:
                    if aminacido.is_disordered() == 0:
                        pass
                    elif aminacido.is_disordered() == 1:
                        alitipo = atomo.get_altloc()
                        alitiparo.append(alitipo)

    # gets just the unique values
    alitiparo = set(alitiparo)
    alitiparo = list(alitiparo)
    try:
        alitiparo.remove(' ')
    except:
        pass

    eligo = output


    # checker for chain type and alternate conformation ID
    class SelectChains(Bio.PDB.Select):
        """ Only accept the specified chains when saving. """

        def __init__(self, chain_letters, atom_altloc_type, model_ids):
            self.chain_letters = chain_letters
            self.atom_altloc_type = atom_altloc_type
            self.model_ids = model_ids
        def accept_model(self, model):
            if (model.get_id() == self.model_ids):
                return (1)
            else:
                return (0)
        def accept_chain(self, chain):
            if (chain.get_id() == self.chain_letters):
                return (1)
            else:
                return (0)
        def accept_atom(self, atom):
            if (atom.get_altloc() == self.atom_altloc_type) \
                    or (atom.get_altloc() == ' '):
                return (1)
            else:
                return (0)

    #writes a file for each chain and alt conf
    outputnames = []
    for modelo in modelaro:
        for cxeno in cxenaro:
            for alitipo in alitiparo:
                cxenonomo = str(cxeno)[10]

                if cxenonomo == " ":
                    novacxenonomo = "X"
                else:
                    novacxenonomo = cxenonomo

                io = PDBIO()
                io.set_structure(strukturo)
                io.save(str(
                            eligo +
                            "_model_" +
                            str(modelo) +
                            "_chain_" +
                            novacxenonomo +
                            "_alt_" +
                            alitipo +
                            ".pdb"
                            ),
                        SelectChains(cxenonomo, alitipo, modelo)
                        )
                # append the names of all the files to a list for easy looping
                # in the next section of code
                outputnames.append(str(
                                        eligo +
                                        "_model_" +
                                        str(modelo) +
                                        "_chain_" +
                                        novacxenonomo +
                                        "_alt_" +
                                        alitipo +
                                        ".pdb"
                                        )
                                    )

    for outputname in outputnames:
        temp = []
        with open(outputname,'r') as endosiero:

            for line in endosiero:
                if (line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM') \
                                          and line[17:20] != 'HOH':
                    new_line = line[0:16] +\
                                 " " + \
                                 line[17:21] + \
                                 "A" + \
                                 line[22:len(line)]
                    eldosiero.write(new_line)
                elif line[17:20] == 'HOH':
                    pass
                else:
                    new_line = line
                    temp.append(new_line)

        os.remove(outputname)
        eldosiero = open(outputname, 'w')
        for line in temp:
            eldosiero.write(line)
        eldosiero.close()

    return outputnames

def chain_order_fixer(outputname):
    #need all these for the loop
    nPrint = False
    caPrint = False
    cPrint = False
    oPrint = False
    caStore = ""
    cStore = ""
    oStore = ""
    nStore = ""
    otherStore = ""
    new_line = ""
    resnum = 1
    prev_resnum = resnum

    temp = []
    with open(outputname,'r') as endosiero:


        temp.append("MODEL        0\n")


        # sets all the chains to be chain A internally, prevents massive
        # errors later
        # removes any indicator of alt_conf within the file
        # now the original chain and the alt conf flag, and the model id only exist
            # in the filename (and thus later in the legend)
        for line in endosiero:
            if line[0:6] == 'ATOM  ':

                if line[13:16] == "N  ":
                    nStore = line[0:16] + \
                             " " + \
                             line[17:21] + \
                             "A" + \
                             line[22:len(line)]
                elif line[13:16] == "CA ":
                    caStore = line[0:16] + \
                              " " + \
                              line[17:21] + \
                              "A" + \
                              line[22:len(line)]
                elif line[13:16] == "C  ":
                    cStore = line[0:16] + \
                             " " + \
                             line[17:21] + \
                             "A" + \
                             line[22:len(line)]
                elif line[13:16] == "O  ":
                    oStore = line[0:16] + \
                             " " + \
                             line[17:21] + \
                             "A" + \
                             line[22:len(line)]
                elif line[13:16] != "N  " \
                        and line[13:16] != "CA " \
                        and line[13:16] != "C  " \
                        and line[13:16] != "O  ":
                    otherStore = otherStore + \
                                 line[0:16] + \
                                 " " + \
                                 line[17:21] + \
                                 "A" + \
                                 line[22:len(line)]

                resnum = int(line[22:26])


                #ensures backbone order
                if nPrint == False \
                        and caPrint == False \
                        and cPrint == False and \
                        oPrint == False:
                    new_line = new_line + nStore
                    nPrint = True
                elif nPrint == True \
                        and caPrint == False \
                        and cPrint == False \
                        and oPrint == False:
                    new_line = new_line + caStore
                    caPrint = True
                elif nPrint == True \
                        and caPrint == True \
                        and cPrint == False \
                        and oPrint == False:
                    new_line = new_line + cStore
                    cPrint = True
                elif nPrint == True \
                        and caPrint == True \
                        and cPrint == True \
                        and oPrint == False:
                    new_line = new_line + oStore
                    oPrint = True

                if resnum != prev_resnum:
                    temp.append(new_line + otherStore)
                    nPrint = False
                    aPrint = False
                    cPrint = False
                    oPrint = False

                    prev_resnum = resnum
                    otherStore = ""
                    new_line = ""

    os.remove(outputname)
    eldosiero = open(outputname, 'w')
    for line in temp:
        eldosiero.write(line)
    eldosiero.write("TER   \n")
    eldosiero.write("ENDMDL")
    eldosiero.close()

# this function checks for chain breaks and gaps in the model
def backbone_scan(pdb):
    # declaring
    resnum = 0
    old_resnum = 0
    atom_lines = False
    atom_order = ''
    filein = open(pdb, "r")
    num_gaps = 0
    for line in filein:
        try:
            if line[0:4] == "ATOM":
                atom_lines = True
                resnum = int(line[22:26])
                if resnum == old_resnum:
                    if line[13:16] == "CA ":
                        atom_order=atom_order+(line[13:16])
                    if line[13:16] == "C  ":
                        atom_order=atom_order+(line[13:16])
                    if line[13:16] == "O  ":
                        atom_order=atom_order+(line[13:16])
                if resnum != old_resnum:
                    if line[13:16] == "N  ":
                        atom_order=atom_order+(line[13:16])
                        if (resnum - old_resnum != 1) and old_resnum != 0:
                            if options.permissive == False \
                                    and options.semipermissive == 0:
                                return True
                            elif options.semipermissive  > 0 \
                                    and options.permissive == False:
                                num_gaps += 1
                        old_resnum = resnum
                    # this would mean it's missing backbone N and thus probably
                    # missing all backbone atoms
                    else:
                        if options.permissive == False \
                                and options.semipermissive == 0:
                            return True
                        elif options.semipermissive  > 0 \
                                and options.permissive == False:
                            num_gaps += 1
        except:
            pass

    #checks to ensure that the only member of this set is one list
    # (should look like this [N,CA,C,O])
    atom_list = atom_order.split()
    all_atoms = [(atom_list[4*x],
                  atom_list[4*x+1],
                  atom_list[4*x+2],
                  atom_list[4*x+3]
                  ) for x in range(len(atom_list)/4)
                  ]
    all_atoms = set(all_atoms)
    all_atoms = list(all_atoms)

    if len(all_atoms) > 1:
        #print pdb, all_atoms
        if options.permissive == False \
                and options.semipermissive == 0:
            return True
        elif options.semipermissive  > 0 and options.permissive == False:
            # if there are more than three types of backbone order
            if len(all_atoms) > (options.semipermissive + 1) \
                or len(all_atoms) > 6:
                return True

    # if semi-permissive is enabled, will remove structures with more than 3
    # gaps
    if num_gaps > options.semipermissive:
        return True

    filein.close()

    # removes files with no atoms
    # this is needed if your input is a file with many models
    # AND alternate chains or conformations
    if atom_lines == False:
        print("No ATOM lines in model " +
              str(pdb)[0:(len(str(pdb))-8)] +
              ", removing from final ensemble."
              )
        return True


# this is a function that will separate all the clusters into their own
# pdb_files
def cluster_sep():
    pdb = open(outputname_all_overlay, 'r')
    atoms = []
    for line in pdb:
        if line[0:4] == "ATOM" \
                or line[0:6] == "HETATM" \
                or line[0:6] == "MODEL " \
                or line[0:6] == "CONECT" \
                or line[0:3] == "END" \
                or line[0:6] == "MASTER" \
                or line[0:3] == "TER" \
                or line[0:6] == "ENDMDL" \
                or line[0:5] == "HELIX" \
                or line[0:5] == "SHEET" \
                or line[0:6] == "REMARK":
            atoms = atoms + [line,]



    # dictionary of all the values to use for coloring
    group_dict = {}
    groups = open('groups_legend.tsv','r')
    # read every line except the first (the header)
    groups_lines = groups.readlines()[1:]
    groups.close()
    for line in groups_lines:
        # if this residue already has a color value stored
        if (int(line.split()[0])) in group_dict:
            pass # do nothing
        else:
            group_dict[int(line.split()[0])] = int(line.split()[1])
            # otherwise store a value in it

    formatted_dict = {}

    for key in group_dict:
        # ensure there are two digits after the decimal place
        replacement = "%0.2f" % (group_dict[key],)
        # ensures the replacement value is the correct number of digits
        while len(replacement) != 5:
            replacement = " " + replacement
        # stores this value for use, should be a five character string
        # like so: 00.10
        formatted_dict[key] = replacement

    # list of new lines that has the b-factors replaced
    new_atoms = []
    key = None

    for line in atoms:

        try:
            key = int(line[11:15])
            new_atoms = new_atoms + [line,]
        except:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                new_atoms = new_atoms + [
                                            str(line[0:60]) +
                                            " " +
                                            str(formatted_dict[key]) +
                                            str(line[67:]),
                                        ]
            else:
                new_atoms = new_atoms + [line,]

    group_list = []
    for line in new_atoms:
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            group_list = group_list + [line[61:66],]
    group_set = set(group_list)
    group_list = list(group_set)

    for key in group_list:


        output = open("group"+ str(int(float(key))) +".pdb", "w")

        old_line = None
        for line in new_atoms:
            try:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    group = line[61:66]
            except:
                pass
            try:
                try:
                    if group == key \
                            and (old_line[0:5] == "MODEL" \
                            or old_line[0:3] == "TER"):
                        output.write(old_line)
                except:
                    pass
                if group == key \
                        and (line[0:5] != "MODEL" \
                        and line[0:3] != "TER"):
                    output.write(line)
            except:
                pass
            old_line = line

        output.close()

# this is a function that will separate all the clusters into their own
# pdb_files
def cluster_sep_non_auto():
    pdb = open(outputname_all_overlay, 'r')
    atoms = []
    for line in pdb:
        if line[0:4] == "ATOM" \
                or line[0:6] == "HETATM" \
                or line[0:6] == "MODEL " \
                or line[0:6] == "CONECT" \
                or line[0:3] == "END" \
                or line[0:6] == "MASTER" \
                or line[0:3] == "TER" \
                or line[0:6] == "ENDMDL" \
                or line[0:5] == "HELIX" \
                or line[0:5] == "SHEET" \
                or line[0:6] == "REMARK":
            atoms = atoms + [line,]

    # dictionary of all the values to use for coloring
    group_dict = {}
    group_dict["excluded"] = " 0.00"
    group_dict["m"] = " 1.00"
    group_dict["n"] = " 2.00"
    group_dict["both"] = " 3.00"

    # list of new lines that has the b-factors replaced
    new_atoms = []
    key = None
    model_num = None
    no_group_N = False
    one_in_both = False

    for line in atoms:
        if line[0:5] == "MODEL":
            # get model number
            model_num = int(line[11:15])

            # if there is a group N
            try:
                if (model_num in options.groupm and \
                    model_num in options.groupn):
                    key = "both"
                    one_in_both = True
                elif model_num in options.groupm:
                    key = "m"
                elif model_num in options.groupn:
                    key = "n"
                else:
                    key = "excluded"
            # if there is only group M
            except:
                no_group_N = True
                if model_num in options.groupm:
                    key = "m"
                else:
                    key = "excluded"

            new_atoms = new_atoms + [line,]
        # if not a MODEL line
        else:
            # if ATOM line replace b factor
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                new_atoms = new_atoms + [
                                            str(line[0:60]) +
                                            " " +
                                            str(group_dict[key]) +
                                            str(line[67:]),
                                        ]
            # otherwise keep the line
            else:
                new_atoms = new_atoms + [line,]
    if no_group_N == False and one_in_both == False:
        group_list = [" 0.00", " 1.00", " 2.00"]
    if no_group_N == True:
        group_list = [" 0.00", " 1.00"]
    if one_in_both == True:
        group_list = [" 0.00", " 1.00", " 2.00", " 3.00"]

    # print a pdb file for each set of models
    for key in group_list:
        if key == " 1.00":
            group_status = "M"
        elif key == " 2.00":
            group_status = "N"
        elif key == " 3.00":
            group_status = "BOTH"
        elif key == " 0.00":
            group_status = "NONE"

        output = open("group" + "_" + group_status +".pdb", "w")

        old_line = None
        for line in new_atoms:
            # read the b factor
            try:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    group = line[61:66]
            except:
                pass
            # if the b factor is the kind we are looking for (ie 2 or 1)
            try:
                try:
                    if group == key \
                            and (old_line[0:5] == "MODEL" \
                            or old_line[0:3] == "TER"):
                        output.write(old_line)
                except:
                    pass
                if group == key \
                        and (line[0:5] != "MODEL" \
                        and line[0:3] != "TER"):
                    output.write(line)
            except:
                pass
            old_line = line

        output.close()
    return(group_list)

# this is a function that will set b factors based on intergroup lodr if
# available, otherwise it will just use group m lodr
def pdb_b_fac(group_list):


    for key in group_list:
        if key == " 1.00":
            group_status = "M"
        elif key == " 2.00":
            group_status = "N"
        elif key == " 3.00":
            group_status = "BOTH"
        elif key == " 0.00":
            group_status = "NONE"



        # open the output and remove it, in order to rewrite our version of it
        new_atoms = []
        with open("group" + "_" + group_status +".pdb", 'r') as pdb:



            value_dict = {}
            # get the intergroup lodr, or the group m lodr
            for resid in eelocal_dict["group_m_lodr"]:
                try:
                    if eelocal_dict\
                               ["group_m_lodr"]\
                               [resid] is not None:
                        value_dict[resid] = float(eelocal_dict\
                                                  ["inter_group_lodr"]\
                                                  [resid]
                                                  )
                except:
                    if eelocal_dict\
                                ["group_m_lodr"]\
                                [resid] is not None:
                        value_dict[resid] = float(eelocal_dict\
                                                  ["group_m_lodr"]\
                                                  [resid]
                                                  )
            # 60 - 66 is the b-factors
            # list of all the ATOM lines
            atoms = []
            # actually it's getting all the lines
            for line in pdb:
                if line[0:4] == "ATOM" or \
                       line[0:6] == "HETATM" or \
                       line[0:6] == "MODEL " or \
                       line[0:6] == "CONECT" or \
                       line[0:3] == "END" or \
                       line[0:6] == "MASTER" or \
                       line[0:3] == "TER" or \
                       line[0:6] == "ENDMDL" or \
                       line[0:5] == "HELIX" or \
                       line[0:5] == "SHEET" or \
                       line[0:6] == "REMARK":
                    atoms = atoms + [line,]

            # figure out what the maximum value in the color values is
            maximum = max(value_dict, key=lambda i: value_dict[i])
            maximum = value_dict[maximum]

            # dictionary of normalized values
            norm_dict = {}

            for key in value_dict:
                # normalize to the maximum
                value = value_dict[key] / maximum
                # ensure there are two digits after the decimal place
                replacement = "%0.2f" % (value,)
                # ensures the replacement value is the correct number of digits
                while len(replacement) != 5:
                    replacement = " " + replacement
                # stores this value for use, should be a five character string like so:
                # " 00.10"
                norm_dict[key] = replacement

            # list of new lines that has the b-factors replaced
            for line in atoms:
                try:
                    key = int(line[22:26])
                    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        if key in norm_dict:
                            # replaces the value
                            new_atoms = new_atoms + [
                                        str(line[0:60]) +
                                            " " +
                                            str(norm_dict[key]) +
                                            str(line[67:]),
                                        ]
                        else:
                            new_atoms = new_atoms + [
                                            str(line[0:60]) +
                                            " 00.00" +
                                            str(line[67:]),
                                        ]
                except:
                    new_atoms = new_atoms + [line,]
        # write the output
        os.remove("group" + "_" + group_status +".pdb")
        output = open("group" + "_" + group_status +".pdb", "w")
        for line in new_atoms:
            output.write(line)
        output.close()

# a funtion to get the distance of any two atoms, given coordinates
def get_dist(atom_y, atom_x):
    xd = atom_y[0] - atom_x[0]
    yd = atom_y[1] - atom_x[1]
    zd = atom_y[2] - atom_x[2]
    distance = math.sqrt(xd*xd + yd*yd + zd*zd)
    return distance






# function to calculate the the LODR score for a specific pair of structures,
# for a specific residue
def eelocal(x,y,resnum):
    # set the superimposer
    sup = SVDSuperimposer()
    # lists of atoms, ref = fixed atoms, alt = moving atoms
    ref_atoms = []
    alt_atoms = []


    # this try statement will fail at any chain breaks as the atoms won't exist
    # if debugging, removing this is a good place to start
        # if you do, tell it to skip the first and last residue manually
    try:
        # get atoms
        ref_atoms.append(structure[x]["A"][resnum - 1]["CA"])
        ref_atoms.append(structure[x]["A"][resnum - 1]["C"])
        ref_atoms.append(structure[x]["A"][resnum - 1]["O"])
        ref_atoms.append(structure[x]["A"][resnum]["N"])
        ref_atoms.append(structure[x]["A"][resnum]["CA"])

        alt_atoms.append(structure[y]["A"][resnum - 1]["CA"])
        alt_atoms.append(structure[y]["A"][resnum - 1]["C"])
        alt_atoms.append(structure[y]["A"][resnum - 1]["O"])
        alt_atoms.append(structure[y]["A"][resnum]["N"])
        alt_atoms.append(structure[y]["A"][resnum]["CA"])


        # get the coords for each atom in the atom lists
        l=len(ref_atoms)
        ref_atoms_coord=np.zeros((l, 3))
        alt_atoms_coord=np.zeros((l, 3))
        for i in range(0, len(ref_atoms)):
            ref_atoms_coord[i]=ref_atoms[i].get_coord()
            alt_atoms_coord[i]=alt_atoms[i].get_coord()


        # superimpose the coords, and then get the rotation matrix and
        # translation
        sup.set(ref_atoms_coord, alt_atoms_coord)
        sup.run()
        rot, tran = sup.get_rotran()


        # lists of atoms to get the LODR scores for
        atoms_x = []
        atoms_y = []

        atoms_x.append(structure[x]["A"][resnum]["C"])
        atoms_x.append(structure[x]["A"][resnum]["O"])
        atoms_x.append(structure[x]["A"][resnum + 1]["N"])
        atoms_x.append(structure[x]["A"][resnum + 1]["CA"])


        atoms_y.append(structure[y]["A"][resnum]["C"])
        atoms_y.append(structure[y]["A"][resnum]["O"])
        atoms_y.append(structure[y]["A"][resnum + 1]["N"])
        atoms_y.append(structure[y]["A"][resnum + 1]["CA"])


        # get coords
        l=len(atoms_x)
        atoms_x_coord=np.zeros((l, 3))
        atoms_y_coord=np.zeros((l, 3))
        for i in range(0, len(atoms_x)):
            atoms_x_coord[i]=atoms_x[i].get_coord()
            atoms_y_coord[i]=atoms_y[i].get_coord()


        # make sure the types are correct, this is straight from the
        # Bio source code
        rot=rot.astype('f')
        tran=tran.astype('f')

        # list of transformed moving atoms
        trans_atom_y_coord = []
        # transform each atom
        for atom in atoms_y_coord:
            trans_atom = np.dot(atom, rot) + tran
            trans_atom_y_coord.append(trans_atom)


        # these four distances
        dist_c = get_dist(trans_atom_y_coord[0],atoms_x_coord[0])
        dist_o = get_dist(trans_atom_y_coord[1],atoms_x_coord[1])
        dist_n = get_dist(trans_atom_y_coord[2],atoms_x_coord[2])
        dist_ca = get_dist(trans_atom_y_coord[3],atoms_x_coord[3])


        #calculate rmsd of the distance for this residue
        dists2 = list()
        dists2.append(dist_c*dist_c)
        dists2.append(dist_o*dist_o)
        dists2.append(dist_n*dist_n)
        dists2.append(dist_ca*dist_ca)
        d2sum = sum(dists2)
        x = d2sum / 4
        lodr = math.sqrt(x)

        #done
        return lodr


    # return Nonetype if any of these atoms cannot be accessed
    except:
        return

# function fo calculate rmsd for a dictionary, will return the same key but
# only one value (the rmsd)
def get_rmsd(atom_dist_dict):

    rmsd_dict = {}

    for key in atom_dist_dict:
        try:
            #calculate rmsd
            dists2 = list()
            for dist in atom_dist_dict[key]:
                dists2.append(dist*dist)
            d2sum = sum(dists2)
            x = d2sum / len(atom_dist_dict[key])
            rmsd = math.sqrt(x)
            rmsd_dict[key] = rmsd

        # this exception here for the lodr scores, which will have some
        # Nonetype values which can't be turned into a RMS
        except:
            rmsd_dict[key] = None

    return rmsd_dict

# same as above basically but for getting the mean
def get_mean(atom_dist_dict):

    mean_dict = {}

    for key in atom_dist_dict:
        try:
            avg = np.mean(atom_dist_dict[key])
            mean_dict[key] = avg
        except:
            mean_dict[key] = None
    return mean_dict


# function to get both the minimum value in a dict[key], and the index from
# that list that gave the min
def get_min(atom_dist_dict):
    # nested dict that has these two things in it
    min_info_dict = NestedDict()
    min_dict = {}
    index_dict = {}
    for entry in atom_dist_dict:
        # this will actually store the minimum and the index of
        # the minimum as tmp
        tmp = min(enumerate(atom_dist_dict[entry]), key=itemgetter(1))
        min_dict[entry] = tmp[1]
        index_dict[entry] = tmp[0]
    min_info_dict["min"] = min_dict
    min_info_dict["index"] = index_dict
    # resturns the nested dict of the min info
    return min_info_dict

# this is the final alignment funtion, which will align all models to the
# first model, using only the common consensus core atoms
def final_aligner(outputname, atoms_to_ignore):
    # first model
    ref_model = structure[0]

    # every other model
    for alt_model in structure:
        all_atom_counter = 0
        counter = 0
        ref_atoms = []
        alt_atoms = []

        # honestly, this looping structure is probably a little
        # obsolete now that I have the chain always set to 'A'
        for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
            for ref_res, alt_res in zip(ref_chain, alt_chain) :
                ref_resid = ref_res.id[1]
                alt_resid = alt_res.id[1]

                for ref_atom, alt_atom in zip(ref_res, alt_res):
                    all_atom_counter += 1
                    ref_atomid = ref_atom.id
                    ref_atom_deets = [ref_resid,ref_atomid]
                    alt_atomid = alt_atom.id
                    alt_atom_deets = [alt_resid,alt_atomid]
                    
                    # need this check because sometimes they are out of order
                    if ref_atom_deets == alt_atom_deets:
                        if not (ref_atom_deets in atoms_to_ignore):
                            counter += 1
                            ref_atoms.append(ref_res[ref_atom.id])
                            alt_atoms.append(alt_res[alt_atom.id]) 

        #Align these paired atom lists:
        super_imposer = Superimposer()

        # exit if all atoms are flagged for removal
        # (ie. if the sets of atoms to superimpose are empty)
        try:
            super_imposer.set_atoms(ref_atoms, alt_atoms)
        except:
            print "\n\n\nCould not superimpose final overlay."
            print("There are no consensus atoms at this cutoff distance that "
                  "are common to all structures!!!\nPlease try again"
                  " with a less strict dcut value.\n\n\n"
                  "ANY eeGLOBAL RESULTS WILL NOT BE ACCURATE UNTIL A "
                  "LESS STRICT DISTANCE CUTOFF IS USED.\n\n\n"
                  )
            return True

        # pretty pointless I suspect
        if ref_model.id == alt_model.id :
            #Check for self/self get zero RMS, zero translation
            #and identity matrix for the rotation.
            assert \
            np.abs(super_imposer.rms) < 0.0000001
            assert \
            np.max(np.abs(super_imposer.rotran[1])) < 0.000001
            assert \
            np.max(np.abs(super_imposer.rotran[0]) - np.identity(3)) < 0.000001
        else :
            #Update the structure by moving all the atoms in
            #this model (not just the ones used for the alignment)
            super_imposer.apply(alt_model.get_atoms())
        percent_core = 100 * (float(counter) / float(all_atom_counter))
        print "RMS(model %i, model %i) = %0.2f (Core contained %0.1f percent of total atoms.)" \
             % (ref_model.id, alt_model.id, super_imposer.rms, percent_core)


    # save the final structure
    io=Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(outputname)

    # fix the pdb file format, otherwise all the models are teated as
    # seperate objects by pymol rather than seperate models
    temp = []
    with open(outputname, 'r') as filein:


        counter = 0
        for line in filein:
            if line[0:3] == 'END' and line[0:6] != 'ENDMDL':
                pass

            # Ensures that the final model numbers in the overlayed file match
            # the legend output by prepare_input
            # this is a critical function
            # otherwise the user will have no way to link the final overlay and
            # statistics back to the original files

            elif line[0:6] == 'MODEL':
                modeltag = str(counter)
                while len(modeltag) < 9:
                    modeltag = " " + modeltag
                if len(modeltag) == 9:
                    modeltag = "MODEL" + modeltag + "\n"
                # write model line
                temp.append(modeltag)
                counter += 1
            else:
                temp.append(line)
    os.remove(outputname)
    fileout = open(outputname, 'w')
    for line in temp:
        fileout.write(line)
    fileout.write("END")
    fileout.close()
    return False

# This is the alignment that will run first when doing the iterative pairwise
# alignments.  This differs from the other aligners in that it will not bother
# to check for core atoms it will just use all atoms to create a pairwise
# alignment.
def first_aligner(x,y):

    # this should be pretty obvious by now
    ref_model = structure[x]
    alt_model = structure[y]

    all_atom_counter = 0
    counter = 0
    ref_atoms = []
    alt_atoms = []
    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res in zip(ref_chain, alt_chain) :
            ref_resid = ref_res.id[1]
            alt_resid = alt_res.id[1]
            # include all atoms (ie. there are no testing criterea here)
            for ref_atom, alt_atom in zip(ref_res, alt_res):
                
                ref_atomid = ref_atom.id
                ref_atom_deets = [ref_resid,ref_atomid]
                alt_atomid = alt_atom.id
                alt_atom_deets = [alt_resid,alt_atomid]
                
                # need this check because sometimes they are out of order
                if ref_atom_deets == alt_atom_deets:
                             
                    ref_atoms.append(ref_res[ref_atom.id])
                    alt_atoms.append(alt_res[alt_atom.id])
                    all_atom_counter += 1
                    counter += 1
                

    #Align these paired atom lists:
    super_imposer = Superimposer()
    # exit if all atoms are flagged for removal (ie. if the sets of atoms
    # to superimpose are empty)
    try:
        super_imposer.set_atoms(ref_atoms, alt_atoms)
    except:
        print( "Initial alignment of model " +
                str(x) +
                " and model " +
                str(y) +
                " failed. Possibly there are no atoms here."
               )
        return

    if ref_model.id == alt_model.id :
        #Check for self/self get zero RMS, zero translation
        #and identity matrix for the rotation.
        assert np.abs(super_imposer.rms) < 0.0000001
        assert np.max(np.abs(super_imposer.rotran[1])) < 0.000001
        assert np.max(np.abs(super_imposer.rotran[0]) - np.identity(3)) \
               < 0.000001
    else :
        #Update the structure by moving all the atoms in
        #this model (not just the ones used for the alignment)
        super_imposer.apply(alt_model.get_atoms())

    # The alignment does not need to return anything, as it's primary function
    # is just to update the coordinates of model y in the pair



# calculates distance per atom for a pair
def per_atom_distance(x,y):
    ref_model = structure[x]
    alt_model = structure[y]
    distance_dict = {}

    for chain in ref_model:
        for res in chain:
            for atom in res:
                # distances can be calculated in BioPython using the subtract
                # function
                distance = structure[x][chain.id][res.id[1]][atom.id] - \
                           structure[y][chain.id][res.id[1]][atom.id]

                # this id_string will be the key, a string with resnum and
                # atomtype connected by a : for easy splitting later
                # this ensures that all equivalent atoms will have the same
                # key, and thus can be compared by the dictionaries

                id_string = str(res.id[1]) + ":" + str(atom.id)
                distance_dict[id_string] = distance

    # return a dictionary of the atom distances
    return distance_dict


# this is the pairwise aligner, which will be called again and again until
# convergence is reached.
def pairwise_realigner(x,y, atoms_to_ignore):

    ref_model = structure[x]
    alt_model = structure[y]
    all_atom_counter = 0
    counter = 0
    ref_atoms = []
    alt_atoms = []

    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res in zip(ref_chain, alt_chain) :
                        
            ref_resid = ref_res.id[1]
            alt_resid = alt_res.id[1]

            for ref_atom, alt_atom in zip(ref_res, alt_res):
                all_atom_counter += 1
                ref_atomid = ref_atom.id
                ref_atom_deets = [ref_resid,ref_atomid]
                alt_atomid = alt_atom.id
                alt_atom_deets = [alt_resid,alt_atomid]
                
                # need this check because sometimes they are out of order
                if ref_atom_deets == alt_atom_deets:
                    if not (ref_atom_deets in atoms_to_ignore):
                        counter += 1
                        ref_atoms.append(ref_res[ref_atom.id])
                        alt_atoms.append(alt_res[alt_atom.id]) 


    #Align these paired atom lists:
    super_imposer = Superimposer()
    # exit if all atoms are flagged for removal (ie. if the sets of atoms to
    # superimpose are empty)
    try:
        super_imposer.set_atoms(ref_atoms, alt_atoms)
    except:
        print("There are no consensus atoms at this cutoff distance that are "
              "common to model %i and %i.\nExiting... Please try again with "
              "a less strict dcut value." % (ref_model.id,
                                             alt_model.id
                                             )
              )
        return

    if ref_model.id == alt_model.id :
        #Check for self/self get zero RMS, zero translation
        #and identity matrix for the rotation.
        assert np.abs(super_imposer.rms) < 0.0000001
        assert np.max(np.abs(super_imposer.rotran[1])) < 0.000001
        assert np.max(np.abs(super_imposer.rotran[0]) - np.identity(3)) \
               < 0.000001
    else :
        #Update the structure by moving all the atoms in
        #this model (not just the ones used for the alignment)
        super_imposer.apply(alt_model.get_atoms())


    # again, does not need to return anything. It is enough to update the
    # coords of model y each time based on the new alignment


def avg_calc_coords(coords1, coords2):
    "Return avg deviations between coords1 and coords2."
    diff = coords1 - coords2
    l = coords1.shape[0]
    return abs((sum(sum(diff)) / l))



# this function calculates the rmsd of the atoms in the core, for any given
# pair of models
# structure is VERY similar to the aligner functions
def get_rms_core(x,y, atoms_to_ignore):
    ref_model = structure[x]
    alt_model = structure[y]
    ref_atoms = []
    alt_atoms = []

    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res in zip(ref_chain, alt_chain) :
            ref_resid = ref_res.id[1]
            alt_resid = alt_res.id[1]

            for ref_atom, alt_atom in zip(ref_res, alt_res):
                ref_atomid = ref_atom.id
                ref_atom_deets = [ref_resid,ref_atomid]
                alt_atomid = alt_atom.id
                alt_atom_deets = [alt_resid,alt_atomid]
                
                if ref_atom_deets == alt_atom_deets:
                    if not (ref_atom_deets in atoms_to_ignore):
                        ref_atoms.append(ref_atom.coord)
                        alt_atoms.append(alt_atom.coord)           

    # actually uses the superimposer SVDsuperimposer which is internal to the
    # Superimposer() function
    sup = SVDSuperimposer()
    # has to deal directly with the coords, nothing else  
    coord1 = np.array(ref_atoms)
    coord2 = np.array(alt_atoms)
    
    # exit if all atoms are flagged for removal (ie. if the sets of atoms to
    # superimpose are empty)
    try:
        if options.avg == False:
            # internal code to calculate rms, not reccomended for users
                # you don't know me
            rms_core = sup._rms(coord1, coord2)
        else:
            rms_core = avg_calc_coords(coord1, coord2)
    except:
        # this elegent error message needs to be reworked when I rewrite all
        # the log write statements.
        # print "RMS GET FAIL"
        return 0.0

    # return the rmsd of the core atoms for this pair
    return round(rms_core,3)


# exactly as above, only it will do so for all atoms
def get_rms_non_core(x,y,atoms_outside_core):

    ref_model = structure[x]
    alt_model = structure[y]
    ref_atoms = []
    alt_atoms = []

    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res in zip(ref_chain, alt_chain) :
            ref_resid = ref_res.id[1]
            alt_resid = alt_res.id[1]

            for ref_atom, alt_atom in zip(ref_res, alt_res):
                ref_atomid = ref_atom.id
                ref_atom_deets = [ref_resid,ref_atomid]
                alt_atomid = alt_atom.id
                alt_atom_deets = [alt_resid,alt_atomid]
                
                if ref_atom_deets == alt_atom_deets:
                    if (ref_atom_deets in atoms_outside_core):
                        ref_atoms.append(ref_atom.coord)
                        alt_atoms.append(alt_atom.coord) 
                        
    # actually uses the superimposer SVDsuperimposer which is internal to the
    # Superimposer() function
    sup = SVDSuperimposer()
    # has to deal directly with the coords, nothing else
    coord1 = np.array(ref_atoms)
    coord2 = np.array(alt_atoms)

    # exit if all atoms are flagged for removal (ie. if the sets of atoms to
    # superimpose are empty)
    try:
        if options.avg == False:
            # internal code to calculate rms, not reccomended for users
                # you don't know me
            rms_non_core = sup._rms(coord1, coord2)
        else:
            rms_non_core = avg_calc_coords(coord1, coord2)
    except:
        # this elegent error message needs to be reworked when I rewrite all
        # the log write statements.
        # print "RMS-NON-CORE GET FAIL"
        return 0.0

    # return the rmsd of the core atoms for this pair
    return round(rms_non_core,3)

# function that for a given pair, will create a list of all atoms that are
# beyond the dcut value away from each other
def dcut_atom_checker(x,y):

    atoms_to_ignore = list()
    all_atoms_counter = 0
    for chain in structure[x]:
        for res in chain:
            for atom in res:
                atom1 = structure[x][chain.id][res.id[1]][atom.id]
                atom2 = structure[y][chain.id][res.id[1]][atom.id]
                distance = atom1 - atom2
                if distance < dcut:
                    pass
                else:
                    atom_deets = [res.id[1],atom.id]
                    atoms_to_ignore.append(atom_deets)
                all_atoms_counter += 1
    
    return atoms_to_ignore, all_atoms_counter

# this function works similarly to eePrep in the prepare_input script.
# it will take the list of atoms to ignore, and get only the unique entries,
# ie, not duplicates. It will return a list of lists that contain [resnum,
# atomtype]
def atom_selector(atoms_to_ignore):

    remove_residues_unformatted = []
    for key in atoms_to_ignore:
        for atom in atoms_to_ignore[key]:
            # it has to be a string that gets appended, in order to be
            # hashable, in order to ensure that only unique values are kept
            remove_residues_unformatted.append(
                                                str(atom[0]) +
                                                "@" +
                                                str(atom[1])
                                                )
    remove_residues_unformatted = set(remove_residues_unformatted)
    remove_residues_unformatted = list(remove_residues_unformatted)

    # now we seperate the string and turn it back into a useful list
    remove_residues = []
    for atom in remove_residues_unformatted:
        resid = int(atom.split("@")[0])
        atomid = atom.split("@")[1]
        remove_residues.append([resid,atomid])

    return remove_residues


def to_single(triple):
    one_letter = {'ALA':'A',
                  'ARG':'R',
                  'ASN':'N',
                  'ASP':'D',
                  'CYS':'C',
                  'GLN':'Q',
                  'GLU':'E',
                  'GLY':'G',
                  'HIS':'H',
                  'ILE':'I',
                  'LEU':'L',
                  'LYS':'K',
                  'MET':'M',
                  'PHE':'F',
                  'PRO':'P',
                  'SER':'S',
                  'THR':'T',
                  'TRP':'W',
                  'TYR':'Y',
                  'VAL':'V',
                 }
    return (one_letter[triple])




################################################################################
#### Prepare Input Files for Analysis ####
################################################################################


def prepare_input(options):

    if not os.path.exists(options.pwd):
        os.makedirs(options.pwd)

    # mark start time

    startTime = time.time()
    # Here the real code begins:


    pdbs = []
    for pdb in options.input:
        pdbs.append(pdb)
    prepped_files_list = []

    # change to the working directory
    os.chdir(options.pwd)

    # runs xray-prep

    print("\n##################\n\nFormatting and seperating..."
          "\n\n##################\n"
          )
    counter = 0
    for pdb in pdbs:
        print("Formatting " +
              str(pdb) +
              " for input, and seperating into " +
              "model, chain, and altconf files."
              )
        xray_outname = "prepped_" + \
            str(ntpath.basename(pdb)[0:(len(ntpath.basename(pdb))-4)])
        prepped_files_list.extend(xrayprep(pdb, xray_outname))
        counter += 1


    # runs backbone_scan
    if options.permissive == False:
        print("\n##################\n\nScanning backbones for gaps and "
              "missing residues...\nTo prevent this, use option"
              " '--permissive'.\n\n##################\n"
              )

    bad_files = []
    for prepped_file in prepped_files_list:
        # if gap or misordered chain, try to fix it
    #    if backbone_scan(prepped_file) == True:
    #        chain_order_fixer(prepped_file)
        # if simple repair didn't work, remove the file
        if backbone_scan(prepped_file) == True:
            print("Removing file: " +
                  str(prepped_file) +
                  " from analysis, due to a gap,\nor atoms being in the " +
                  "wrong order in the file (this breaks the ensemblator).\n"
                  )
            bad_files.append(prepped_file)




    good_files = set(prepped_files_list) - set(bad_files)
    good_files = list(good_files)


    # run this before eePrep in order to ensure that the files are good for prep
    if options.align == True:

        from Bio.PDB import *
        from Bio.Align.Applications import MuscleCommandline
        from Bio import SeqIO

        io = PDBIO()
        pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)

        align_counter = 1
        while align_counter != 0:
            # get all the sequences for the pdb files
            seqs = open("sequences.fasta", "w")
            for pdb in good_files:
                structure = pdb_reader.get_structure("temp", pdb)

                seq = ""
                for residue in structure[0]["A"]:
                    try:
                        seq = seq + str(to_single(residue.get_resname()))
                    except:
                        # it's probably DNA or RNA rather than protein
                        if residue.get_resname() in ['A','C','T','G']:
                            seq = seq + str(residue.get_resname())
                        else:
                            seq = seq + "X"
                        
                seqs.write(">" + str(pdb) + "\n" + str(seq) + "\n")
            seqs.close()

            # command line args for MUSCLE
            cline = MuscleCommandline(input="sequences.fasta",
                                      out="muscle_align.fasta"
                                      )
            print( "\n##################\n\n" + \
                      "Running MUSCLE alignment program with the following " +\
                      "command: \n" + str(cline) +\
                      "\n\n##################\n"
                  )

            # RUN MUCSLE
            stdout, stderr = cline()
            os.remove("sequences.fasta")

            aligned_dict = {}

            aligned_file = open("muscle_align.fasta", "rU")
            # fill a dictionary with just the filenames and the aligned sequences
            for element in SeqIO.parse(aligned_file, "fasta"):
                aligned_dict[element.id] = str(element.seq)
            aligned_file.close()


            # ensure that the seq id is close enough to use
            # if a model is defined

            template_name = "prepped_" + \
                            options.template + \
                            "_model_" + \
                            str(options.model) + \
                            "_chain_" + \
                            str(options.chain) + \
                            "_alt_A.pdb"


            #empty the good_files list, to refill with the newly good files
            # ie the files with a certain percent id
            good_files = []
            align_counter = 0
            for key in aligned_dict:
                try:
                    template = aligned_dict[template_name]
                except:
                    print("Couldn't load template file, probably it was removed"
                             " after backbone scanning. Try running with the"
                             " --permissive option on, or pick a different "
                             "template. Alternativly it is not in the set"
                             " of files you are building in the ensemble (this"
                             " is required), or the filename does not end in "
                             "'.pdb'.\n\n\n Alternativly, you have not"
                             " selected a correct file, chain, or model ID."
                             " If there is no chain ID, use 'X'. These chain"
                             " IDs are case sensitive, so use A rather than a."
                            )
                tester = aligned_dict[key]
                #calculate percent id
                matches = 0.0
                for pos in range(0,(len(template))):
                    if tester[pos] == template[pos]:
                        matches += 1.0
                percent_id = matches/float(len(template))
                percent = options.percent
                # if the user entered a fraction or a percent, use a fraction
                if percent > 1.0:
                    percent = percent / 100.0
                # keep the file if it's good
                if percent_id >= percent:
                    good_files.append(key)
                else:
                    message = "Removing file: " + \
                              key + \
                              " due to having less than " + \
                              str(percent * 100) + "%" +\
                              " identity to the template."
                    print message
                    # we removed files, so redo the alignmnet
                    align_counter += 1

        # for each file
        for key in good_files:
            structure = pdb_reader.get_structure("temp", key)
            pos_counter = 0
            index_counter = -1
            old_pos = 0
            old_index = 0
            # for each residue in the structure
            for residue in structure[0]["A"]:
                old_pos = pos_counter
                old_index = index_counter
                pos_counter += 1
                index_counter += 1
                # as long as reading in a gap, increase the residue number
                try:
                    while aligned_dict[key][index_counter] == '-':
                        pos_counter += 1
                        index_counter += 1
                # will break if the last residue in the sequence is a gap
                except:
                    pass
                # if not a gap in the alignment
                # set the residue id to be the universal position
                residue.id = (' ', pos_counter, ' ')
            # save the newly numbered structure
            io.set_structure(structure)
            io.save(key)






    print("\n##################\n\nCombining all files into an ensemble ready " +
          "for use with ensemblator by ensuring that all atoms match in all " +
          "structures...\n\n##################\n"
          )


    # runs eeprep
    legend_dict = eeprep(good_files, bad_files)
    legend_file = open('model_legend.tsv', 'w')
    legend_file.write('model\tfilename\n')
    counter = 0
    for key in legend_dict:
        if legend_dict[key] == "REMOVED FROM FINAL ENSEMBLE":
            legend_file.write("NA" + "\t" + str(legend_dict[key]) + '\n')
        else:
            legend_file.write(str(counter) +
                              "\t" +
                              str(legend_dict[key])[8:len(str(legend_dict[key]))] +
                              '\n')
            counter += 1
    legend_file.close()

    #cleaning up
    for prepped_file in prepped_files_list:
        try:
            os.remove(prepped_file)
        except:
            pass
    print("\n##################\n\nWrote file: " +
          str(options.output) +
          "\n\n##################\n"
          )


    # align all the structures to make the output look nicer

    print "\n##################\n\nCalculating overlay...\n\n##################\n"
    try:
        aligner(options.output)
    except:
        err_str=("The overlay failed. This usually means that there were no"
              " atoms in the ensemble, which in turn usually means that afer"
              " removing the atoms not common to all structures, you were left"
              " with nothing. \n\n\n Sometimes this occurs because all of"
              " the structures you input have gaps, and you are not"
              " allowing structures with gaps. Scroll up through"
              " the log to see if that is the case, and if so, try allowing"
              " some or all gaps.\n\n\n"
              )
        print(err_str)



    # time stuff
    endTime = time.time()
    workTime =  endTime - startTime
    print "\n\nCompleted in: " + str(workTime) + " seconds.\n"




################################################################################
#### Analyze ####
################################################################################

def analyze(options):


    if not os.path.exists(options.pwd):
        os.makedirs(options.pwd)


    # mark start time
    startTime = time.time()

    # read in the input
    pdb = options.input

    # change to the working directory
    os.chdir(options.pwd)

    global dcut
    dcut = options.dcut

    pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)
    outputname = "global_overlay_" + str(dcut) + ".pdb"
    # use this in the auto analysis, cluster_sep - this is a lazy solution
    global outputname_all_overlay
    outputname_all_overlay = "global_overlay_" + str(dcut) + ".pdb"

    # this structure variable is constantly referenced by functions.
    # Do not mess around with it lightly.
    global structure
    structure = pdb_reader.get_structure("temp", pdb)

    print("Iterativly aligning pairwise structures until convergence of cutoff "
          "accepted atoms is reached..."
          )


    # first output file, contains info about all the pairs. Used for clustering
    # in the --auto analysis
    pairwise_file = open("pairwise_analysis.tsv", 'w')
    if options.avg == False:
        pairwise_file.write("model_X\tmodel_Y\tcore_percent\trmsd_non_core\trmsd_core\tdis_score\n")
    else:
        pairwise_file.write("model_X\tmodel_Y\tcore_percent\tavg_dev_non_core\tavg_dev_core\tdis_score\n")

    # take each pairwise alignment and realign until dcut satisfied atoms converge

    # get the list of models to generate the pairs
    model_list = list()
    for model in structure:
        model_list.append(model.id)

    atoms_to_ignore = {}
    # list of all pairs
    pairwise_list = combinations(model_list, r = 2)
    # iterate accross each pair in order, ignoring duplicates and self pairs
    # (they aren't in the pairwise list)
    for x,y in pairwise_list:
        # just any random starting condition. I don't think this is really needed
        # but hey.
        atoms = 1
        atoms2 = 500
        counter = 0

        # do first alignment of the pair, using all atoms
        first_aligner(x,y)
        # will realign over and over until the list of kept atoms are exactly
        # identical before and after overlay
        while atoms != atoms2:
            counter += 1
            # do first dcut check
            atoms, all_atoms = dcut_atom_checker(x,y)
            # here the key in atoms_to_ignore is a unique string for each pair,
            atoms_to_ignore[str(x) + "," + str(y)] = atoms
            if atoms == 0:
                # this would mean that the cutoff is too strict
                print "No atoms converge at cutoff."
                break
            # re-align with only core atoms, giving the pairwise aligner the same
            # list of atoms to ignore that is specific to this pair
            pairwise_realigner(x,y, atoms_to_ignore[str(x) + "," + str(y)])
            # now the atoms that pass the check are here in atoms2, and if they
            # aren't all identical to atoms, then the loop will repeat
            atoms2, all_atoms = dcut_atom_checker(x,y)


        # now that a convergent core is found, calculate the stats for this pair
        rms_core = get_rms_core(x,y,atoms_to_ignore[str(x) + "," + str(y)])
        rms_non_core = get_rms_non_core(x,y, atoms_to_ignore[str(x) + "," + str(y)])


        core_percent = round((float(all_atoms) - float(len(atoms2))) / float(all_atoms), 3)
        dis_score = round(math.pow(rms_core, core_percent) * \
                    math.pow(rms_non_core, (1 - core_percent)), 3)
        
        
        # output information to table, tab separated
        pairwise_file.write(
                            str(x) +
                            "\t" +
                            str(y) +
                            "\t" +
                            str(core_percent)+
                            "\t" +
                            str(rms_non_core) +
                            "\t" +
                            str(rms_core) +
                            "\t" +
                            str(dis_score) +
                            "\n"
                            )


    # self-self pairs, needed to get a complete array later
    # all these stats should be zero
    for modelid in model_list:
        # include all self-self pairs in the output table to make matrix
        # generation more complete if one desires
        pairwise_file.write(
                            str(modelid) +
                            "\t" +
                            str(modelid) +
                            "\t" +
                            str(0)+
                            "\t" +
                            str(0)+
                            "\t" +
                            str(0) +
                            "\t" +
                            str(0) +
                            "\n"
                            )

    # generate common-core aligned file
    print("Identifying common convergant atoms in all pairwise structures, and "
          "realigning original ensemble using only common cutoff accepted atoms..."
          )
    # get the common elements from the dict of all atoms to ignore from all
    # structures
    removal_list = atom_selector(atoms_to_ignore)
    counter = 0
    for atom in removal_list:
        counter += 1
    print(str(counter) +
          " non-consensus atoms removed from final aligner process."
          )

    # now re-align the original structure, against the first model, using only the
    # COMMON core atoms atoms
    # will be true if there are no common core atoms
    error_check = final_aligner(outputname, removal_list)

    if error_check != True:
        print("Wrote final overlay as: " + outputname)


    print("Saved number of rejected atoms per pair, RMSD of all atoms per pair,"
          "and RMSD of only core atoms per pair, as 'pairwise_analysis.tsv'."
          )
    pairwise_file.close()

    # this is all the eeGlobal and eeLocal stuff that will run if the m and n
    # group are set manually
    if options.auto == False:
        # NOW DOING eeGLOBAL stuff
        # calulate RMS level of intra-ensemble variation for each atom among
        # m structures
        print "Calculating Intra Group M RMSD:"
        pairwise_list = combinations(options.groupm, r = 2)
        all_dist_dict = {}

        # for each pair of models in group m
        for x,y in pairwise_list:
            # get the distances per atom
            distance_dict = per_atom_distance(x,y)
            for key in distance_dict:
                # append the value to a list of pairwise distances, by atom
                if key in all_dist_dict:
                    all_dist_dict[key].append(distance_dict[key])
                else:
                    all_dist_dict[key] = list()
                    all_dist_dict[key].append(distance_dict[key])

        # get either the avg or the rmsd of these values
        if options.avg == False:
            group_m_rmsd = get_rmsd(all_dist_dict)
        else:
            group_m_rmsd = get_mean(all_dist_dict)

        # try statement here is so that it will except if there is only a group M
        # defined
        # calulate RMS level of intra-ensemble variation for each atom among n
        # structures
        try:
            # same as group m
            pairwise_list = combinations(options.groupn, r = 2)
            all_dist_dict = {}
            print "Calculating Intra Group N RMSD:"
            for x,y in pairwise_list:
                distance_dict = per_atom_distance(x,y)
                for key in distance_dict:
                    if key in all_dist_dict:
                        all_dist_dict[key].append(distance_dict[key])
                    else:
                        all_dist_dict[key] = list()
                        all_dist_dict[key].append(distance_dict[key])
            if options.avg == False:
                group_n_rmsd = get_rmsd(all_dist_dict)
            else:
                group_n_rmsd = get_mean(all_dist_dict)

            # calulate RMS level of inter-ensemble variation for each atom
            # between m and n structures
            # same as group m and n, but will be done on the pairs that are m-n
            print "Calculating Inter Group RMSD:"
            inter_list = [x for x in itertools.product(
                                                        options.groupm,
                                                        options.groupn
                                                        )
                          ]

            all_dist_dict = {}
            for x,y in inter_list:
                distance_dict = per_atom_distance(x,y)
                for key in distance_dict:
                    if key in all_dist_dict:
                        pass
                    else:
                        all_dist_dict[key] = list()
                    all_dist_dict[key].append(distance_dict[key])
            if options.avg == False:
                inter_rmsd = get_rmsd(all_dist_dict)
            else:
                inter_rmsd = get_mean(all_dist_dict)

            print "Calculating Closest Approach Distance:"
            # now get the closest approach info
            closest_approach_info = get_min(all_dist_dict)
            closest_approach_index =  closest_approach_info["index"]
            closest_approach = closest_approach_info["min"]

        except:
            pass

        ### combine them all into a well formatted dictionary for ease of access
        eeglobal_dict = NestedDict()
        for key in group_m_rmsd:
            # remember the key format is "resid:atomtype"
            resid = int(key.split(":")[0])
            atomid = key.split(":")[1]
            eeglobal_dict\
                           ["group_m_rmsd"]\
                           [resid]\
                           [atomid]\
                           = group_m_rmsd[key]
            try:
                eeglobal_dict\
                               ["group_n_rmsd"]\
                               [resid]\
                               [atomid]\
                               = group_n_rmsd[key]
                eeglobal_dict\
                               ["inter_group_rmsd"]\
                               [resid]\
                               [atomid]\
                               = inter_rmsd[key]
                eeglobal_dict\
                               ["closest_approach"]\
                               [resid]\
                               [atomid]\
                               = closest_approach[key]
                eeglobal_dict\
                               ["closest_approach_index"]\
                               [resid]\
                               [atomid]\
                               = str(inter_list[closest_approach_index[key]])
            except:
                pass

        # sort them for nicer output
        # this ensures that the lists/tables are in a nicer order for humans
        # this is used later for looping, so it's needed
        resid_list = []
        for key in eeglobal_dict["group_m_rmsd"]:
            resid_list.append(key)
        resid_list = set(resid_list)
        resid_list = list(sorted(resid_list))

        atomid_list = []
        for resid in eeglobal_dict["group_m_rmsd"]:
            for key in eeglobal_dict["group_m_rmsd"][resid]:
                atomid_list.append(key)
        atomid_list = set(atomid_list)
        atomid_list = list(sorted(atomid_list))



        #### eeLocal calculations
        # all methods here are as above, only they only use resnum as keys,
        # so are less complicated
        # they also calculate lodr rather than rmsd, etc


        print "Calculating LODR scores for group M:"
        pairwise_list = combinations(options.groupm, r = 2)
        all_lodr_dict = {}
        lodr_dict = {}

        for x,y in pairwise_list:
            for resnum in resid_list:
                lodr_dict[resnum] = eelocal(x,y, int(resnum))

                if resnum in all_lodr_dict:
                    all_lodr_dict[resnum].append(lodr_dict[resnum])
                else:
                    all_lodr_dict[resnum] = list()
                    all_lodr_dict[resnum].append(lodr_dict[resnum])

        if options.avg == False:
            group_m_lodr = get_rmsd(all_lodr_dict)
        else:
            group_m_lodr = get_mean(all_lodr_dict)

        ### calulate LODR among n structures
        try:
            pairwise_list = combinations(options.groupn, r = 2)
            print "Calculating LODR scores for group N:"
            all_lodr_dict = {}
            lodr_dict = {}

            for x,y in pairwise_list:
                for resnum in resid_list:
                    lodr_dict[resnum] = eelocal(x,y, int(resnum))

                    if resnum in all_lodr_dict:
                        all_lodr_dict[resnum].append(lodr_dict[resnum])
                    else:
                        all_lodr_dict[resnum] = list()
                        all_lodr_dict[resnum].append(lodr_dict[resnum])

            if options.avg == False:
                group_n_lodr = get_rmsd(all_lodr_dict)
            else:
                group_n_lodr = get_mean(all_lodr_dict)
            ### calulate LODR between m and n structures
            print "Calculating Inter Group LODR:"

            inter_list = [x for x in itertools.product(
                                                        options.groupm,
                                                        options.groupn
                                                        )
                          ]

            all_lodr_dict = {}
            lodr_dict = {}
            for x,y in inter_list:
                for resnum in resid_list:
                    lodr_dict[resnum] = eelocal(x,y, int(resnum))

                    if resnum in all_lodr_dict:
                        all_lodr_dict[resnum].append(lodr_dict[resnum])
                    else:
                        all_lodr_dict[resnum] = list()
                        all_lodr_dict[resnum].append(lodr_dict[resnum])
            if options.avg == False:
                inter_group_lodr = get_rmsd(all_lodr_dict)
            else:
                inter_group_lodr = get_mean(all_lodr_dict)

            print "Calculating Minimum LODR between M and N at each residue:"

            minimum_lodr_info = get_min(all_lodr_dict)
            minimum_lodr_index =  minimum_lodr_info["index"]
            minimum_lodr = minimum_lodr_info["min"]

        except:
            pass


        # same as eeGlobal
        eelocal_dict = NestedDict()
        for resid in resid_list:
            try:
                eelocal_dict\
                              ["group_m_lodr"]\
                              [resid]\
                              = group_m_lodr[resid]
                try:
                    eelocal_dict\
                                  ["group_n_lodr"]\
                                  [resid]\
                                  = group_n_lodr[resid]
                    eelocal_dict\
                                  ["inter_group_lodr"]\
                                  [resid]\
                                  = inter_group_lodr[resid]
                    eelocal_dict\
                                  ["minimum_lodr"]\
                                  [resid]\
                                  = minimum_lodr[resid]
                    eelocal_dict\
                                  ["minimum_lodr_index"]\
                                  [resid]\
                                  = str(inter_list[minimum_lodr_index[resid]])
                except:
                    pass
            except:
                pass


        # print all the output tables
        print "Organizing output, and saving files:"

        # output tables
        eeglobal_out = open("eeGlobal_out.tsv", 'w')

        # header
        eeglobal_out.write("res_id" +
                           "\t" +
                           "atom_id" +
                           "\t" +
                           "intra_m_rmsd" +
                           "\t" +
                           "intra_n_rmsd" +
                           "\t" +
                           "inter_group_rmsd" +
                           "\t" +
                           "closest_approach" +
                           "\t" +
                           "closest_approach_pair" +
                           "\t" +
                           "included_in_core" +
                           "\n"
                           )
        for resid in resid_list:
            for atomid in atomid_list:
                # if the atom exists at this residue
                if atomid in eeglobal_dict["group_m_rmsd"][resid]:
                    # if there is another group, otherwise print NA in
                    # those columns
                    if "group_n_rmsd" in eeglobal_dict:
                        # if the atom is in the core or was removed
                        if [int(resid), atomid] in removal_list:
                            eeglobal_out.write(str(resid) +
                                               "\t" +
                                               str(atomid) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["group_m_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["group_n_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["inter_group_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["closest_approach"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["closest_approach_index"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               "False" +
                                               "\n"
                                               )
                        else:
                            eeglobal_out.write(str(resid) +
                                               "\t" +
                                               str(atomid) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["group_m_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["group_n_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["inter_group_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["closest_approach"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["closest_approach_index"]
                                                   [resid]
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               "True" +
                                               "\n"
                                               )
                    else:
                        if [int(resid), atomid] in removal_list:
                            eeglobal_out.write(str(resid) +
                                               "\t" +
                                               str(atomid) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["group_m_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "False" +
                                               "\n"
                                               )
                        else:
                            eeglobal_out.write(str(resid) +
                                               "\t" +
                                               str(atomid) +
                                               "\t" +
                                               str(eeglobal_dict\
                                                   ["group_m_rmsd"]\
                                                   [resid]\
                                                   [atomid]
                                                   ) +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "NA" +
                                               "\t" +
                                               "True" +
                                               "\n"
                                               )
                else:
                    pass
        eeglobal_out.close()
        print "Output saved in file 'eeGlobal_out.tsv'."


        eelocal_out = open("eeLocal_out.tsv", 'w')
        eelocal_out.write("res_id"
                          "\t"
                          "intra_m_rms_lodr"
                          "\t"
                          "intra_n_rms_lodr"
                          "\t"
                          "inter_group_rms_lodr"
                          "\t"
                          "minimum_lodr"
                          "\t"
                          "mimimum_lodr_pair"
                          "\n"
                          )
        for resid in resid_list:
            if "group_n_lodr" in eelocal_dict:
                eelocal_out.write(str(resid) +
                                  "\t" +
                                  str(eelocal_dict["group_m_lodr"][resid]) +
                                  "\t" +
                                  str(eelocal_dict["group_n_lodr"][resid]) +
                                  "\t" +
                                  str(eelocal_dict\
                                      ["inter_group_lodr"]\
                                      [resid]
                                      ) +
                                  "\t" +
                                  str(eelocal_dict["minimum_lodr"][resid]) +
                                  "\t" +
                                  str(eelocal_dict\
                                      ["minimum_lodr_index"]\
                                      [resid]
                                      ) +
                                  "\n"
                                  )
            else:
                eelocal_out.write(str(resid) +
                                  "\t" +
                                  str(eelocal_dict["group_m_lodr"][resid]) +
                                  "\t" +
                                  "NA" +
                                  "\t" +
                                  "NA" +
                                  "\t" +
                                  "NA" +
                                  "\t" +
                                  "NA" +
                                  "\n"
                                  )
        eelocal_out.close()
        print "Output saved in file 'eeLocal_out.tsv'."

        ### PLOTTING ###
        # don't know what this does
        minorLocator = MultipleLocator(1)
        # colors
        red = '#e41a1c'
        blue = '#377eb8'
        green = '#4daf4a'
        purple = '#decbe4'
        yellow = '#ffff33'

        print "Plotting eeGlobal:"
        ## eeGlobal plot
        # we are jsut interested in showing the average backbone rmsd here
        backbone_intra_m_rmsd = {}
        for resid in range(min(resid_list),(max(resid_list)+1)):
            rmsds = []
            for atomid in atomid_list:
                if atomid == "N" or \
                           atomid == "CA" or \
                           atomid == "C" or \
                           atomid == "O":
                    rmsds.append(eeglobal_dict["group_m_rmsd"][resid][atomid])
            try:
                backbone_intra_m_rmsd[resid] = np.mean(rmsds)
            except:
                backbone_intra_m_rmsd[resid] = None

        try:
            backbone_intra_n_rmsd = {}
            for resid in range(min(resid_list),(max(resid_list)+1)):
                rmsds = []
                for atomid in atomid_list:
                    if atomid == "N" or \
                               atomid == "CA" or \
                               atomid == "C" or \
                               atomid == "O":
                        rmsds.append(eeglobal_dict["group_n_rmsd"][resid][atomid])
                try:
                    backbone_intra_n_rmsd[resid] = np.mean(rmsds)
                except:
                    backbone_intra_n_rmsd[resid] = None

            backbone_inter_rmsd = {}
            for resid in range(min(resid_list),(max(resid_list)+1)):
                rmsds = []
                for atomid in atomid_list:
                    if atomid == "N" or \
                               atomid == "CA" or \
                               atomid == "C" or \
                               atomid == "O":
                        rmsds.append(eeglobal_dict\
                                     ["inter_group_rmsd"]\
                                     [resid]\
                                     [atomid]
                                     )
                try:
                    backbone_inter_rmsd[resid] = np.mean(rmsds)
                except:
                    backbone_inter_rms[resid] = None
            backbone_closest = {}
            for resid in range(min(resid_list),(max(resid_list)+1)):
                rmsds = []
                for atomid in atomid_list:
                    if atomid == "N" or \
                               atomid == "CA" or \
                               atomid == "C" or \
                               atomid == "O":
                        rmsds.append(eeglobal_dict\
                                     ["closest_approach"]\
                                     [resid]\
                                     [atomid]
                                     )
                try:
                    backbone_closest[resid] = np.mean(rmsds)
                except:
                    backbone_closest[resid] = None
        except:
            pass


        if options.avg == False:
            title = "eeGLOBAL_dcut=" + str(dcut)
            plt.figure()
            plt.plot(backbone_intra_m_rmsd.keys(),
                       backbone_intra_m_rmsd.values(),
                       blue,
                       label="Group M RMSD",
                       linewidth=1.5
                       )
            try:
                plt.plot(backbone_intra_n_rmsd.keys(),
                           backbone_intra_n_rmsd.values(),
                           green,
                           label="Group N RMSD",
                           linewidth=1.5
                           )
                plt.plot(backbone_inter_rmsd.keys(),
                           backbone_inter_rmsd.values(),
                           purple,
                           label="Between groups RMSD",
                           linewidth=1.5
                           )
                plt.plot(backbone_closest.keys(),
                           backbone_closest.values(),
                           red,
                           label="Closest Approach",
                           linewidth=1.5
                           )
            except:
                pass
            plt.xlabel("Residue Number")
            plt.ylabel("Backbone rmsd")
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                         loc=3,
                         ncol=2,
                         mode="expand",
                         borderaxespad=0.
                         )
            ax = plt.gca()
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            fig = plt.gcf()
            fig.canvas.set_window_title(title)
            plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
            plt.close()
            print "eeGlobal plot saved as '" + title + ".png'."


            print "Plotting eeLocal:"
            ## eeLocal plot

            title = "eeLocal"
            plt.figure()
            plt.plot(eelocal_dict["group_m_lodr"].keys(),
                       eelocal_dict["group_m_lodr"].values(),
                       blue,
                       label="Group M RMS-LODR",
                       linewidth=1.5
                       )
            try:
                plt.plot(eelocal_dict["group_n_lodr"].keys(),
                           eelocal_dict["group_n_lodr"].values(),
                           green,
                           label="Group N RMS-LODR",
                           linewidth=1.5
                           )
                plt.plot(eelocal_dict["inter_group_lodr"].keys(),
                           eelocal_dict["inter_group_lodr"].values(),
                           purple,
                           label="Inter-group RMS-LODR",
                           linewidth=1.5
                           )
                plt.plot(eelocal_dict["minimum_lodr"].keys(),
                           eelocal_dict["minimum_lodr"].values(),
                           red,
                           label="Minimum LODR",
                           linewidth=1.5
                           )
            except:
                pass
            plt.xlabel("Residue Number")
            plt.ylabel("RMS-LODR")
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                         loc=3,
                         ncol=2,
                         mode="expand",
                         borderaxespad=0.
                         )
            ax = plt.gca()
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            fig = plt.gcf()
            fig.canvas.set_window_title(title)
            plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
            plt.close()
            print "eeLocal plot saved as '" + title + ".png'."

        # calculated using averages here, so plots should reflect that
        else:

            title = "eeGLOBAL_dcut=" + str(dcut)
            plt.figure()
            plt.plot(backbone_intra_m_rmsd.keys(),
                       backbone_intra_m_rmsd.values(),
                       blue,
                       label="Group M Mean",
                       linewidth=1.5
                       )
            try:
                plt.plot(backbone_intra_n_rmsd.keys(),
                           backbone_intra_n_rmsd.values(),
                           green,
                           label="Group N Mean",
                           linewidth=1.5
                           )
                plt.plot(backbone_inter_rmsd.keys(),
                           backbone_inter_rmsd.values(),
                           purple,
                           label="Between groups Mean",
                           linewidth=1.5
                           )
                plt.plot(backbone_closest.keys(),
                           backbone_closest.values(),
                           red,
                           label="Closest Approach",
                           linewidth=1.5
                           )
            except:
                pass
            #plt.title(title)
            plt.xlabel("Residue Number")
            plt.ylabel("Backbone rmsd")
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                         loc=3,
                         ncol=2,
                         mode="expand",
                         borderaxespad=0.
                         )
            #plt.subplots_adjust(right=.78)
            ax = plt.gca()
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            fig = plt.gcf()
            fig.canvas.set_window_title(title)
            plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
            plt.close()
            print "eeGlobal plot saved as '" + title + ".png'."


            print "Plotting eeLocal:"
            ## eeLocal plot

            title = "eeLocal"
            plt.figure()
            plt.plot(eelocal_dict["group_m_lodr"].keys(),
                       eelocal_dict["group_m_lodr"].values(),
                       blue,
                       label="Group M Mean-LODR",
                       linewidth=1.5
                       )
            try:
                plt.plot(eelocal_dict["group_n_lodr"].keys(),
                           eelocal_dict["group_n_lodr"].values(),
                           green,
                           label="Group N Mean-LODR",
                           linewidth=1.5
                           )
                plt.plot(eelocal_dict["inter_group_lodr"].keys(),
                           eelocal_dict["inter_group_lodr"].values(),
                           purple,
                           label="Inter-group Mean-LODR",
                           linewidth=1.5
                           )
                plt.plot(eelocal_dict["minimum_lodr"].keys(),
                           eelocal_dict["minimum_lodr"].values(),
                           red,
                           label="Minimum LODR",
                           linewidth=1.5
                           )
            except:
                pass
            plt.xlabel("Residue Number")
            plt.ylabel("RMS-LODR")
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                         loc=3,
                         ncol=2,
                         mode="expand",
                         borderaxespad=0.
                         )
            ax = plt.gca()
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            fig = plt.gcf()
            fig.canvas.set_window_title(title)
            plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
            plt.close()
            print "eeLocal plot saved as '" + title + ".png'."

        # this will set the b-factors in the final overlay to be the
        # value of the inter_group b factor. Can be nice for visualizing
        if options.color == True:
            print("Setting b-factors to relative LODR score in final overlay file, "
                  "for use with 'spectrum b' command in pymol."
                  )
            group_list = cluster_sep_non_auto()
            pdb_b_fac(group_list)
        else:
            print "Creating pdb files for each group..."
            cluster_sep_non_auto()

    ###############################################################################
    # AUTOMATIC HERE:
    if options.auto == True:

        # this normalizes the datasets
        # it divides every observation by the stdev of all the observations
        # (this is different than the default whiten function)
        def whiten(obs):
            std_dev = np.std(obs)
            if std_dev != 0.0:
                return obs / std_dev
            else:
                return obs

        # read in the pairwise data file, we gonna build a matrix
        pairwise_file = open("pairwise_analysis.tsv", 'r')

        # read every line except the first (the header)
        lines = pairwise_file.readlines()[1:]
        pairwise_file.close()

        # a nested dictionary for each of these three stats we are interested in
        dis_score = NestedDict()


        # build pairwise dictionaries with all the values of interest
        for line in lines:
            x = int(line.split()[0])
            y = int(line.split()[1])

            dis = float(line.split()[5])
            dis_score[x][y] = dis

        # this is used to generate ranges, it's so I know what the highest number
        # of models is (ie. how many columns and rows I need in my matrix)
        max_y = int(max(list(dis_score.keys())))


        # making a list that will be formatted in such a way that I can turn it
        # into an array or a matrix
        dis_score_array_list = []

        # go from 0 to the max number of models
        for x in range(0,max_y + 1):
            this_x_dis_score = []
            # now do the same for y, now we are going over every pair
            for y in range(0,max_y + 1):

                # fill out the missing values in the array, by duplicating the
                # correct data
                if dis_score[x][y] == {}:
                    dis_score[x][y] = dis_score[y][x]

                # append these dictionary values to the list
                this_x_dis_score.append(dis_score[x][y])

            # now append them all, what the list looks like for each x (ie. row)
            # is this [0,23,43,23,53,654,23] where index 0 is x0y0, index 1 is
            # x0y1 etc.
            dis_score_array_list.append(this_x_dis_score)


        #### generate co-ocur matrix

        co_ocur = NestedDict()

        # go from 0 to the max number of models
        for x in range(0,max_y + 1):
            # now do the same for y, now we are going over every pair
            for y in range(0,max_y + 1):
                co_ocur[x][y] = 0


        print "\nPerforming iterative Affinity Propagation clustering:\n"
        
        # affinity propagation using dis score
        X = np.array(dis_score_array_list)
        # convert from a distance matrix to a similarity matrix
        X = X * -1

        # start with minimum pref
        pref = np.min(X[np.nonzero(X)])


        af = AffinityPropagation(preference = pref,
                                 affinity = "precomputed",
                                 max_iter=2000).fit(X)

        cluster_centers_indices = af.cluster_centers_indices_
        labels = af.labels_
        pref = pref * 1.01

        # now loop without recording until there are less than n clusters
        while len(np.unique(labels)) == len(labels):
            af = AffinityPropagation(preference = pref,
                                     affinity = "precomputed",
                                     max_iter=2000).fit(X)

            cluster_centers_indices = af.cluster_centers_indices_
            labels = af.labels_
            pref = pref * 1.01


        # now record and loop until one cluster
        while len(np.unique(labels)) != 1:
            
            # update co-ocur matrix
            # go from 0 to the max number of models
            for x in range(0,max_y + 1):
                # now do the same for y, now we are going over every pair
                for y in range(0,max_y + 1):
                    
                    if labels[x] == labels[y]:
                        co_ocur[x][y] = co_ocur[x][y] + 1
            

            af = AffinityPropagation(preference = pref,
                                     affinity = "precomputed",
                                     max_iter=2000).fit(X)

            cluster_centers_indices = af.cluster_centers_indices_
            labels = af.labels_
            pref = pref * 1.01

            
            

        ################## k-means now

        print "\nPerforming iterative K-means clustering:\n"

        # get whitened matrixes of these bad boys!
        dis_score_asmatrix = whiten(np.asmatrix(dis_score_array_list))

        # loop from 2 to n-minus-one/2 clusters, some number i times each
        for k in range(2, len(labels)):
            
            i = 0
            
            while i < 10:
            
            
                codebook, dis_score_distortion = kmeans(dis_score_asmatrix, k)
                labels, dist = vq(dis_score_asmatrix, codebook)
                dis_score_best_distortion = dis_score_distortion
                
                i += 1
                
                # update co-ocur matrix
                # go from 0 to the max number of models
                for x in range(0,max_y + 1):
                    # now do the same for y, now we are going over every pair
                    for y in range(0,max_y + 1):
                        
                        if labels[x] == labels[y]:
                            co_ocur[x][y] = co_ocur[x][y] + 1


        print "\nPerforming Agglomerative clustering of ensemble of clusters:\n"


        # now cluster the co_ocur matrix
        co_ocur_array_list = []
        # go from 0 to the max number of models
        for x in range(0,max_y + 1):
            this_x_co_ocur = []
            # now do the same for y, now we are going over every pair
            for y in range(0,max_y + 1):
                # append these dictionary values to the list
                this_x_co_ocur.append(co_ocur[x][y])

            # now append them all, what the list looks like for each x (ie. row)
            # is this [0,23,43,23,53,654,23] where index 0 is x0y0, index 1 is
            # x0y1 etc.
            co_ocur_array_list.append(this_x_co_ocur)

        X = np.array(co_ocur_array_list)

        # get the max number of clusters to search for from the user
        # optional. Default is 6 (as declared in the option.parser at the top)
        max_clust = options.maxclust + 1

        # same as above
        sil_score = -1.0
        sil_score_best = -1.0
        labels_best = []

        for k in range(2,max_clust):
            
            labels = AgglomerativeClustering(k, affinity = "euclidean", linkage = "complete").fit_predict(X)
            sil_score = metrics.silhouette_score(X, labels, metric='euclidean')
            
            # now check again and set the value
            if sil_score > sil_score_best or sil_score_best == -1.0:
                labels_best = labels
                sil_score_best = sil_score


        num_clust = len(np.unique(labels))
                
        sil_scores = metrics.silhouette_samples(X, labels, metric='euclidean')
        sil_scores_out = open("sil_scores.tsv", "w")
        counter = 0
        sil_scores_out.write("id" + "\t" + 
                             "cluster" + "\t" + 
                             "sil_score" + "\n")
        for label in labels:
            sil_scores_out.write(str(counter) + "\t" + 
                                 str(label) + "\t" + 
                                 str(sil_scores[counter]) + "\n")
            counter += 1


        print "\nThere are " + str(num_clust) + " clusters, with a mean " + \
              "silhouette score of " + str(sil_score) + "." 
        print "Sillhouette Scores saved in 'sil_scores.tsv'\n"

        best_code = labels        



        # get the pairwise list of all the clusters compared against each other,
        # as we did with x,y values above
        
        cluster_combos = combinations(range(int(num_clust)),2)

        # do analysis for each group vs each group!!
        # we need a legend in order to know which models belong to which group
        # used to include this in the names of the plots, but they got too long

        groups_legend = open("groups_legend.tsv", "w")

        for groupm,groupn in cluster_combos:

            group_m = []
            group_n = []
            counter = 0
            for model in best_code:
                if model == groupm:
                    group_m.append(counter)
                    groups_legend.write(str(counter) + "\t" + str(groupm) + "\n")
                if model == groupn:
                    group_n.append(counter)
                    groups_legend.write(str(counter) + "\t" + str(groupn) + "\n")
                counter += 1
            
            
            # redo but reverse groups if group m only has one member,
            # this prevents empty tables and graphs
            if len(group_m) == 1:
                
                group_m = []
                group_n = []
                counter = 0
                for model in best_code:
                    if model == groupn:
                        group_m.append(counter)
                        groups_legend.write(str(counter) + "\t" + str(groupn) + "\n")
                    if model == groupm:
                        group_n.append(counter)
                        groups_legend.write(str(counter) + "\t" + str(groupm) + "\n")
                    counter += 1
                outputname = "_Group_" + str(groupn) + "_Group_" + str(groupm)
            else:
                outputname = "_Group_" + str(groupm) + "_Group_" + str(groupn)
            
            
            print "\nComparing clusters " + str(groupm) + " and " + str(groupn) + ".\n"

            # NOW DOING eeGLOBAL stuff
            # same as above, in the non-auto section.
            # for detailed comments, see that section of the code

            # calulate RMS level of intra-ensemble variation for each atom
            # among m structures
            try:
                print "Calculating Intra Group M RMSD:"
                pairwise_list = combinations(group_m, r = 2)
                all_dist_dict = {}

                for x,y in pairwise_list:
                    distance_dict = per_atom_distance(x,y)
                    for key in distance_dict:
                        if key in all_dist_dict:
                            all_dist_dict[key].append(distance_dict[key])
                        else:
                            all_dist_dict[key] = list()
                            all_dist_dict[key].append(distance_dict[key])

                if options.avg == False:
                    group_m_rmsd = get_rmsd(all_dist_dict)
                else:
                    group_m_rmsd = get_mean(all_dist_dict)
            except:
                pass

            # calulate RMS level of intra-ensemble variation for each atom
            # among n structures
            try:

                pairwise_list = combinations(group_n, r = 2)
                all_dist_dict = {}
                print "Calculating Intra Group N RMSD:"
                for x,y in pairwise_list:
                    distance_dict = per_atom_distance(x,y)
                    for key in distance_dict:
                        if key in all_dist_dict:
                            all_dist_dict[key].append(distance_dict[key])
                        else:
                            all_dist_dict[key] = list()
                            all_dist_dict[key].append(distance_dict[key])
                if options.avg == False:
                    group_n_rmsd = get_rmsd(all_dist_dict)
                else:
                    group_n_rmsd = get_mean(all_dist_dict)
            except:
                pass
            try:
                # calulate RMS level of inter-ensemble variation for each atom
                # between m and n structures
                print "Calculating Inter Group RMSD:"
                inter_list = [x for x in itertools.product(group_m, group_n)]
                all_dist_dict = {}
                for x,y in inter_list:
                    distance_dict = per_atom_distance(x,y)
                    for key in distance_dict:
                        if key in all_dist_dict:
                            pass
                        else:
                            all_dist_dict[key] = list()
                        all_dist_dict[key].append(distance_dict[key])
                if options.avg == False:
                    inter_rmsd = get_rmsd(all_dist_dict)
                else:
                    inter_rmsd = get_mean(all_dist_dict)
                print "Calculating Closest Approach Distance:"

                closest_approach_info = get_min(all_dist_dict)
                closest_approach_index =  closest_approach_info["index"]
                closest_approach = closest_approach_info["min"]

            except:
                pass

            # combine them all into a well formatted dictionary for ease of
            # access

            eeglobal_dict = NestedDict()


            for key in group_m_rmsd:
                resid = int(key.split(":")[0])
                atomid = key.split(":")[1]
                try:
                    eeglobal_dict\
                                   ["group_m_rmsd"]\
                                   [resid]\
                                   [atomid]\
                                   = group_m_rmsd[key]
                except:
                    pass
                try:
                    eeglobal_dict\
                                   ["group_n_rmsd"]\
                                   [resid]\
                                   [atomid]\
                                   = group_n_rmsd[key]
                except:
                    pass
                try:
                    eeglobal_dict\
                                   ["inter_group_rmsd"]\
                                   [resid]\
                                   [atomid]\
                                   = inter_rmsd[key]
                    eeglobal_dict\
                                   ["closest_approach"]\
                                   [resid]\
                                   [atomid]\
                                   = closest_approach[key]
                    eeglobal_dict\
                                   ["closest_approach_index"]\
                                   [resid]\
                                   [atomid]\
                                   = str(inter_list[closest_approach_index[key]])
                except:
                    pass

            # sort them for nicer output
            resid_list = []
            for key in eeglobal_dict["group_m_rmsd"]:
                resid_list.append(key)
            resid_list = set(resid_list)
            resid_list = list(sorted(resid_list))

            atomid_list = []
            for resid in eeglobal_dict["group_m_rmsd"]:
                for key in eeglobal_dict["group_m_rmsd"][resid]:
                    atomid_list.append(key)
            atomid_list = set(atomid_list)
            atomid_list = list(sorted(atomid_list))



            # eeLocal calculations
            try:
                print "Calculating LODR scores for group M:"
                pairwise_list = combinations(group_m, r = 2)
                all_lodr_dict = {}
                lodr_dict = {}

                for x,y in pairwise_list:
                    for resnum in resid_list:
                        lodr_dict[resnum] = eelocal(x,y, int(resnum))

                        if resnum in all_lodr_dict:
                            all_lodr_dict[resnum].append(lodr_dict[resnum])
                        else:
                            all_lodr_dict[resnum] = list()
                            all_lodr_dict[resnum].append(lodr_dict[resnum])

                if options.avg == False:
                    group_m_lodr = get_rmsd(all_lodr_dict)
                else:
                    group_m_lodr = get_mean(all_lodr_dict)
            except:
                pass
            # calulate LODR among n structures
            try:
                pairwise_list = combinations(group_n, r = 2)
                print "Calculating LODR scores for group N:"
                all_lodr_dict = {}
                lodr_dict = {}

                for x,y in pairwise_list:
                    for resnum in resid_list:
                        lodr_dict[resnum] = eelocal(x,y, int(resnum))

                        if resnum in all_lodr_dict:
                            all_lodr_dict[resnum].append(lodr_dict[resnum])
                        else:
                            all_lodr_dict[resnum] = list()
                            all_lodr_dict[resnum].append(lodr_dict[resnum])

                if options.avg == False:
                    group_n_lodr = get_rmsd(all_lodr_dict)
                else:
                    group_n_lodr = get_mean(all_lodr_dict)
            except:
                pass
            try:
                # calulate LODR between m and n structures
                print "Calculating Inter Group LODR:"
                inter_list = [x for x in itertools.product(group_m, group_n)]

                all_lodr_dict = {}
                lodr_dict = {}
                for x,y in inter_list:
                    for resnum in resid_list:
                        lodr_dict[resnum] = eelocal(x,y, int(resnum))

                        if resnum in all_lodr_dict:
                            all_lodr_dict[resnum].append(lodr_dict[resnum])
                        else:
                            all_lodr_dict[resnum] = list()
                            all_lodr_dict[resnum].append(lodr_dict[resnum])
                if options.avg == False:
                    inter_group_lodr = get_rmsd(all_lodr_dict)
                else:
                    inter_group_lodr = get_mean(all_lodr_dict)

                print("Calculating Minimum LODR between " +
                      "M and N at each residue:"
                      )

                minimum_lodr_info = get_min(all_lodr_dict)
                minimum_lodr_index =  minimum_lodr_info["index"]
                minimum_lodr = minimum_lodr_info["min"]

            except:
                pass



            eelocal_dict = NestedDict()
            for resid in resid_list:
                try:
                    try:
                        eelocal_dict\
                                      ["group_m_lodr"]\
                                      [resid]\
                                      = group_m_lodr[resid]
                    except:
                        pass
                    try:
                        eelocal_dict\
                                      ["group_n_lodr"]\
                                      [resid]\
                                      = group_n_lodr[resid]
                    except:
                        pass
                    try:
                        eelocal_dict\
                                      ["inter_group_lodr"]\
                                      [resid]\
                                      = inter_group_lodr[resid]
                        eelocal_dict\
                                      ["minimum_lodr"]\
                                      [resid]\
                                      = minimum_lodr[resid]
                        eelocal_dict\
                                      ["minimum_lodr_index"]\
                                      [resid]\
                                      = str(
                                          inter_list[minimum_lodr_index[resid]]
                                          )
                    except:
                        pass
                except:
                    pass



            print "Organizing output, and saving files:"

            # output tables
            eeglobal_out = open("eeGlobal_out" + outputname + ".tsv", 'w')
            eeglobal_out.write("res_id"
                               "\t"
                               "atom_id"
                               "\t"
                               "intra_m_rmsd"
                               "\t"
                               "intra_n_rmsd"
                               "\t"
                               "inter_group_rmsd"
                               "\t"
                               "closest_approach"
                               "\t"
                               "closest_approach_pair"
                               "\t"
                               "included_in_core"
                               "\n"
                               )
            for resid in resid_list:
                for atomid in atomid_list:
                    if atomid in eeglobal_dict["group_m_rmsd"][resid]:
                        if "group_n_rmsd" in eeglobal_dict:
                            if [int(resid), atomid] in removal_list:
                                eeglobal_out.write(str(resid) +
                                                   "\t" +
                                                   str(atomid) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["group_m_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["group_n_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["inter_group_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["closest_approach"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["closest_approach_index"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   "False" +
                                                   "\n"
                                                   )
                            else:
                                eeglobal_out.write(str(resid) +
                                                   "\t" +
                                                   str(atomid) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["group_m_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["group_n_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["inter_group_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["closest_approach"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["closest_approach_index"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   "True" +
                                                   "\n"
                                                   )
                        else:
                            if [int(resid), atomid] in removal_list:
                                eeglobal_out.write(str(resid) +
                                                   "\t" +
                                                   str(atomid) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["group_m_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "False" +
                                                   "\n"
                                                   )
                            else:
                                eeglobal_out.write(str(resid) +
                                                   "\t" +
                                                   str(atomid) +
                                                   "\t" +
                                                   str(eeglobal_dict\
                                                       ["group_m_rmsd"]\
                                                       [resid]\
                                                       [atomid]
                                                       ) +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "NA" +
                                                   "\t" +
                                                   "True" +
                                                   "\n"
                                                   )
                    else:
                        pass
            eeglobal_out.close()
            print "Output saved in file 'eeGlobal_out" + outputname + ".tsv'."

            eelocal_out = open("eeLocal_out" + outputname + ".tsv", 'w')
            eelocal_out.write("res_id"
                              "\t"
                              "intra_m_rms_lodr"
                              "\t"
                              "intra_n_rms_lodr"
                              "\t"
                              "inter_group_rms_lodr"
                              "\t"
                              "minimum_lodr"
                              "\t"
                              "mimimum_lodr_pair"
                              "\n"
                              )

            for resid in resid_list:
                if "group_n_lodr" in eelocal_dict:
                    eelocal_out.write(str(resid) +
                                      "\t" +
                                      str(eelocal_dict\
                                          ["group_m_lodr"]\
                                          [resid]
                                          ) +
                                      "\t" +
                                      str(eelocal_dict\
                                          ["group_n_lodr"]\
                                          [resid]
                                          ) +
                                      "\t" +
                                      str(eelocal_dict\
                                          ["inter_group_lodr"]\
                                          [resid]
                                          ) +
                                      "\t" +
                                      str(eelocal_dict\
                                          ["minimum_lodr"]\
                                          [resid]
                                          ) +
                                      "\t" +
                                      str(eelocal_dict\
                                          ["minimum_lodr_index"]\
                                          [resid]
                                          ) +
                                      "\n"
                                      )
                else:
                    eelocal_out.write(str(resid) +
                                      "\t" +
                                      str(eelocal_dict\
                                          ["group_m_lodr"]\
                                          [resid]) +
                                      "\t" +
                                      "NA" +
                                      "\t" +
                                      "NA" +
                                      "\t" +
                                      "NA" +
                                      "\t" +
                                      "NA" +
                                      "\n"
                                      )
            eelocal_out.close()
            print "Output saved in file 'eeLocal_out" + outputname + ".tsv'."

            ### PLOTTING ###

            minorLocator = MultipleLocator(1)
            red = '#e41a1c'
            blue = '#377eb8'
            green = '#4daf4a'
            purple = '#decbe4'
            yellow = '#ffff33'

            print "Plotting eeGlobal:"
            ## eeGlobal plot
            backbone_intra_m_rmsd = {}
            try:
                for resid in range(min(resid_list),(max(resid_list)+1)):
                    rmsds = []
                    for atomid in atomid_list:
                        if atomid == "N" or \
                                   atomid == "CA" or \
                                   atomid == "C" or \
                                   atomid == "O":
                            rmsds.append(eeglobal_dict["group_m_rmsd"][resid][atomid])
                    try:
                        backbone_intra_m_rmsd[resid] = np.mean(rmsds)
                    except:
                        backbone_intra_m_rmsd[resid] = None
            except:
                pass
            try:
                backbone_intra_n_rmsd = {}
                for resid in range(min(resid_list),(max(resid_list)+1)):
                    rmsds = []
                    for atomid in atomid_list:
                        if atomid == "N" or \
                                   atomid == "CA" or \
                                   atomid == "C" or \
                                   atomid == "O":
                            rmsds.append(eeglobal_dict\
                                         ["group_n_rmsd"]\
                                         [resid]\
                                         [atomid]
                                         )
                    try:
                        backbone_intra_n_rmsd[resid] = np.mean(rmsds)
                    except:
                        backbone_intra_n_rmsd[resid] = None
            except:
                pass
            try:
                backbone_inter_rmsd = {}
                for resid in range(min(resid_list),(max(resid_list)+1)):
                    rmsds = []
                    for atomid in atomid_list:
                        if atomid == "N" or \
                                   atomid == "CA" or \
                                   atomid == "C" or \
                                   atomid == "O":
                            rmsds.append(eeglobal_dict\
                                         ["inter_group_rmsd"]\
                                         [resid]\
                                         [atomid]
                                         )
                    try:
                        backbone_inter_rmsd[resid] = np.mean(rmsds)
                    except:
                        backbone_inter_rmsd[resid] = None
                backbone_closest = {}
                for resid in range(min(resid_list),(max(resid_list)+1)):
                    rmsds = []
                    for atomid in atomid_list:
                        if atomid == "N" or \
                                   atomid == "CA" or \
                                   atomid == "C" or \
                                   atomid == "O":
                            rmsds.append(eeglobal_dict\
                                         ["closest_approach"]\
                                         [resid]\
                                         [atomid]
                                         )
                    try:
                        backbone_closest[resid] = np.mean(rmsds)
                    except:
                        backbone_closest[resid] = None
            except:
                pass


            if options.avg == False:
                title = "eeGLOBAL_dcut=" + str(dcut) + outputname
                plt.figure()
                try:
                    plt.plot(backbone_intra_m_rmsd.keys(),
                               backbone_intra_m_rmsd.values(),
                               blue,
                               label="Group M RMSD",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(backbone_intra_n_rmsd.keys(),
                               backbone_intra_n_rmsd.values(),
                               green,
                               label="Group N RMSD",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(backbone_inter_rmsd.keys(),
                               backbone_inter_rmsd.values(),
                               purple,
                               label="Between groups RMSD",
                               linewidth=1.5
                               )
                    plt.plot(backbone_closest.keys(),
                               backbone_closest.values(),
                               red,
                               label="Closest Approach",
                               linewidth=1.5
                               )
                except:
                    pass
                plt.xlabel("Residue Number")
                plt.ylabel("Backbone rmsd")
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                             loc=3,
                             ncol=2,
                             mode="expand",
                             borderaxespad=0.
                             )
                ax = plt.gca()
                ax.xaxis.set_minor_locator(minorLocator)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                fig = plt.gcf()
                fig.canvas.set_window_title(title)
                plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
                plt.close()
                print "eeGlobal plot saved as '" + title + ".png'."


                print "Plotting eeLocal:"
                ## eeLocal plot

                title = "eeLocal" + outputname
                plt.figure()
                try:
                    plt.plot(eelocal_dict["group_m_lodr"].keys(),
                               eelocal_dict["group_m_lodr"].values(),
                               blue,
                               label="Group M RMS-LODR",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(eelocal_dict["group_n_lodr"].keys(),
                               eelocal_dict["group_n_lodr"].values(),
                               green, label="Group N RMS-LODR",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(eelocal_dict["inter_group_lodr"].keys(),
                               eelocal_dict["inter_group_lodr"].values(),
                               purple,
                               label="Inter-group RMS-LODR",
                               linewidth=1.5
                               )
                    plt.plot(eelocal_dict["minimum_lodr"].keys(),
                               eelocal_dict["minimum_lodr"].values(),
                               red,
                               label="Minimum LODR",
                               linewidth=1.5
                               )
                except:
                    pass
                plt.xlabel("Residue Number")
                plt.ylabel("RMS-LODR")
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                             loc=3,
                             ncol=2,
                             mode="expand",
                             borderaxespad=0.
                             )
                ax = plt.gca()
                ax.xaxis.set_minor_locator(minorLocator)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                fig = plt.gcf()
                fig.canvas.set_window_title(title)
                plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
                plt.close()
                print "eeLocal plot saved as '" + title + ".png'."

            else:

                title = "eeGLOBAL_dcut=" + str(dcut) + outputname
                plt.figure()
                try:
                    plt.plot(backbone_intra_m_rmsd.keys(),
                               backbone_intra_m_rmsd.values(),
                               blue,
                               label="Group M Mean",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(backbone_intra_n_rmsd.keys(),
                               backbone_intra_n_rmsd.values(),
                               green,
                               label="Group N Mean",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(backbone_inter_rmsd.keys(),
                               backbone_inter_rmsd.values(),
                               purple,
                               label="Between groups Mean",
                               linewidth=1.5
                               )
                    plt.plot(backbone_closest.keys(),
                               backbone_closest.values(),
                               red,
                               label="Closest Approach",
                               linewidth=1.5
                               )
                except:
                    pass
                plt.xlabel("Residue Number")
                plt.ylabel("Backbone rmsd")
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                             loc=3,
                             ncol=2,
                             mode="expand",
                             borderaxespad=0.
                             )
                ax = plt.gca()
                ax.xaxis.set_minor_locator(minorLocator)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                fig = plt.gcf()
                fig.canvas.set_window_title(title)
                plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
                plt.close()
                print "eeGlobal plot saved as '" + title + ".png'."


                print "Plotting eeLocal:"
                ## eeLocal plot

                title = "eeLocal" + outputname
                plt.figure()
                try:
                    plt.plot(eelocal_dict["group_m_lodr"].keys(),
                               eelocal_dict["group_m_lodr"].values(),
                               blue, label="Group M Mean-LODR",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(eelocal_dict["group_n_lodr"].keys(),
                               eelocal_dict["group_n_lodr"].values(),
                               green,
                               label="Group N Mean-LODR",
                               linewidth=1.5
                               )
                except:
                    pass
                try:
                    plt.plot(eelocal_dict["inter_group_lodr"].keys(),
                               eelocal_dict["inter_group_lodr"].values(),
                               purple,
                               label="Inter-group Mean-LODR",
                               linewidth=1.5
                               )
                    plt.plot(eelocal_dict["minimum_lodr"].keys(),
                               eelocal_dict["minimum_lodr"].values(),
                               red,
                               label="Minimum LODR",
                               linewidth=1.5
                               )
                except:
                    pass
                plt.xlabel("Residue Number")
                plt.ylabel("RMS-LODR")
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                             loc=3,
                             ncol=2,
                             mode="expand",
                             borderaxespad=0.
                             )
                ax = plt.gca()
                ax.xaxis.set_minor_locator(minorLocator)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                fig = plt.gcf()
                fig.canvas.set_window_title(title)
                plt.savefig(title + ".png" , dpi = 400, bbox_inches='tight')
                plt.close()
                print "eeLocal plot saved as '" + title + ".png'."



        # now need to sort the groups legend so that it is human readable
        groups_legend = open("groups_legend.tsv", "r")
        lines = [line for line in groups_legend if line.strip()]
        lines = set(lines)
        lines = list(lines)
        groups_legend.close()
        os.remove("groups_legend.tsv")

        lines_fixed = []
        for line in lines:
            new_line = line.replace('\t', '.')
            new_line = float(new_line)
            lines_fixed.append(new_line)
        lines_fixed.sort()

        groups_legend_sorted = open("groups_legend.tsv", "w")
        groups_legend_sorted.write("model\tgroup\n")
        for line in lines_fixed:
            new_line = str(line)
            new_line = new_line.replace('.','\t')
            groups_legend_sorted.write(new_line + "\n")
        groups_legend_sorted.close()


        # leave this here, it needs groups_legend.tsv to still exist
        print "Creating pdb files for each group..."
        cluster_sep()


        # now if the model_legend.tsv file exists, append the groups to that
        # and remove the group_legend file. Otherwise don't

        #if it exists
        try:
            # handles for the two legends
            mLegend = open("model_legend.tsv", "r")
            gLegend = open("groups_legend.tsv", "r")

            # create lists, from all the lines, one at a time,
            # so that I can iterate over one, but use the same index
            # to refer to the same model from the other
            mLines = []
            gLines = []
            nLines = []
            mismatch = False

            for line in mLegend:
                mLines.append(line)
            for line in gLegend:
                gLines.append(line)


            counter = 0
            # rewrite the legend
            for line in mLines:
                try:
                    # redundent check to ensure the models are the same (checks #)
                    if gLines[counter].split("\t")[0] == line.split("\t")[0]:
                        nLines.append(line.strip() +
                                     "\t" +
                                     gLines[counter].split("\t")[1]
                                     )
                    else:
                        print "Group legend saved seperatly " + \
                              "due to model # mismatch."
                        print "Saved group identity in 'groups_legend.tsv'"
                        mismatch = True
                        break
                    counter += 1
                except:
                    # model and group legend don't contain the same number of
                    # members. Don't mess with either file, just leave it.
                    # user error
                    print "Group legend saved seperatly due to model # mismatch."
                    print "Saved group identity in 'groups_legend.tsv'"
                    mismatch = True
                    break

            # if everything went well, write the new legend. Otherwise leave
            # everything the same
            if mismatch == False:
                os.remove("groups_legend.tsv")
                os.remove("model_legend.tsv")
                output = open("model_legend.tsv", "w")
                for line in nLines:
                    output.write(line)
                output.close()
                print "Saved group identity in 'model_legend.tsv'"

            mLegend.close()
            gLegend.close()

        # no model legend file
        except:
            pass


    # time stuff
    endTime = time.time()
    workTime =  endTime - startTime
    print "\nEverything completed in: " + str(workTime) + " seconds."
    print "\nResults (images and ensembles) saved in " + options.pwd + "\n"




    # Done!




#### GUI STUFF ####

class Utility:
    def __init__(self):
        self.applocation = os.path.dirname(sys.argv[0]);
        self.last_dir_accessed = os.path.realpath(os.environ['HOME']);
        self.dir_to_save = ""

    def select_output_dir(self,label_to_update):
        dir_to_save = tkFileDialog.askdirectory(initialdir = self.last_dir_accessed, title = "Select output directory.");
        dir_to_save_label = "..." + dir_to_save[len(dir_to_save) - 22:]

        self.last_dir_accessed = dir_to_save

        label_to_update["text"] = dir_to_save_label
        return(dir_to_save);

    def move_files(self, from_to_dict, base_path):
        for from_file in from_to_dict:
            if(os.path.isfile(os.path.join(base_path,from_file))):
                shutil.move(os.path.join(base_path,from_file), os.path.join(base_path,from_to_dict[from_file]));

    def select_input_files(self, min_number, label_to_update):
        chosenfiles = tkFileDialog.askopenfilename(initialdir = self.last_dir_accessed, multiple = True, title="Select input pdb files.");
        toRet = list();

        firstname = os.path.basename(chosenfiles[0]);
        lastname = os.path.basename(chosenfiles[len(chosenfiles)-1]);
        label_to_update["text"] = firstname + " ... " + lastname
        if self.dir_to_save != "":
            self.last_dir_accessed = os.path.dirname(self.dir_to_save)
        else:
            self.last_dir_accessed = os.path.dirname(chosenfiles[0])
        toRet = chosenfiles;

        return(toRet);

    def select_input_file(self, label_to_update):
        chosenfile = tkFileDialog.askopenfilename(initialdir = self.last_dir_accessed,
                                                  multiple = False,
                                                  title="Select template file. Dissimilar chains will be removed from the ensemble.");
        toRet = "";
        if(chosenfile != ""):
            firstname = os.path.basename(chosenfile);
            label_to_update["text"] = firstname
            print self.dir_to_save
            if self.dir_to_save != "":
                self.last_dir_accessed = (self.dir_to_save)
            else:
                self.last_dir_accessed = os.path.dirname(chosenfile)
            toRet = chosenfile;
        return(toRet);

    def select_input_ensemble(self, label_to_update):
        chosenfile = tkFileDialog.askopenfilename(initialdir = self.last_dir_accessed,
                                                  multiple = False,
                                                  title="Select input ensemble (must have been prepared using the Ensemblator.");
        toRet = "";
        if(chosenfile != ""):
            firstname = os.path.basename(chosenfile);
            label_to_update["text"] = firstname
            if self.dir_to_save != "":
                self.last_dir_accessed = (self.dir_to_save)
            else:
                self.last_dir_accessed = os.path.dirname(chosenfile)
            toRet = chosenfile;
        return(toRet);


    def prepare_input_command_run(self, options, label_to_update, window_to_update):

        try:
            prepare_input(options)
        except Exception as e:
            tkMessageBox.showerror("Oops.", str(e));

        label_to_update["text"] = "Done!";
        window_to_update.update();

    def analyze_command_run(self, options, label_to_update, window_to_update):

        try:
            analyze(options)
        except Exception as e:
            tkMessageBox.showerror("Oops.", str(e));

        label_to_update["text"] = "Done!";
        window_to_update.update();


class MainRoot:
    def __init__(self):
        self.rootWindow = Tk();
        self.rootWindow.wm_title("The Ensemblator");
        self.setup_gui();
        self.utility = Utility();
        self.rootWindow.mainloop();

    def setup_gui(self):
        self.add_row("Prepare Input",
                     "Prepare ensemble for analysis.",
                     self.prepare,
                     0
                     )
        self.add_row("Analyze",
                     "Analyze prepared ensemble.",
                     self.analyze,
                     1
                     )
        self.add_row("Exit",
                     "Exit application.",
                     self.exit,
                     2
                     )

        t1 = Text(self.rootWindow)
        sys.stdout = TextRedirector(t1, "stdout")
        sys.stderr = TextRedirector(t1, "stderr")
        t1.grid(row = 0, rowspan = 3, column = 2, sticky=E)


    def add_row(self, button_name, label_text, func_to_call, rownum):

        frame1 = Frame(self.rootWindow);
        #frame1.grid(row=0, column = 0);
        fr1_button = Button(self.rootWindow,
                            text = button_name,
                            command = func_to_call
                            )
        fr1_button.grid(row = rownum, column = 0, sticky = E+W);
        fr1_label = Label(self.rootWindow, text = label_text);
        fr1_label.grid(row = rownum, column = 1, sticky = W);

    def prepare(self):
        subwindow = Prepare(self.utility);
    def analyze(self):
        subwindow = Analyze(self.utility);
    def exit(self):
        sys.exit("Have a nice day!")



class Prepare:
    def __init__(self, utility):
        self.utility = utility;
        self.rootWindow = Tk();
        self.rootWindow.wm_title("Prepare Input");
        self.output = StringVar(self.rootWindow);
        self.output.set("prepared_ensemble.pdb");
        self.break_num = IntVar(self.rootWindow);
        self.break_num.set(0);
        self.gap_setting = StringVar(self.rootWindow);
        self.gap_setting.set("None");
        self.align = IntVar(self.rootWindow)
        self.align.set(0)
        self.template = ""
        self.chain = StringVar(self.rootWindow)
        self.chain.set("X")
        self.model = StringVar(self.rootWindow)
        self.model.set("0")
        self.percent = StringVar(self.rootWindow)
        self.percent.set("70")


        self.inputfiles = list();
        self.dir_to_save = "";
        self.setup_gui();
        self.rootWindow.mainloop();

    def setup_gui(self):

        select_button = Button(self.rootWindow,
                               text = "Select Input Files",
                               command = self.select_files
                               )
        select_button.grid(row=0,
                           column = 0,
                           columnspan = 2,
                           sticky = E+W
                           );

        files_selected_label = Label(self.rootWindow,
                                     text = "Files Selected: "
                                     )
        files_selected_label.grid(row = 1,
                                  column = 0,
                                  sticky = W
                                  )
        self.files_selected = Label(self.rootWindow,
                                    text = "None           "
                                    )
        self.files_selected.grid(row = 1,
                                 column = 1,
                                 columnspan = 3,
                                 sticky = W
                                 )


        pwd_button = Button(self.rootWindow,
                            text = "Select Working Directory",
                            command = self.select_pwd
                            )
        pwd_button.grid(row=2,
                        column = 0,
                        columnspan = 2,
                        sticky = E+W
                        )

        pwd_selected_label = Label(self.rootWindow,
                                   text = "Working Directory: "
                                   )
        pwd_selected_label.grid(row = 3,
                                column = 0,
                                sticky = W
                                )
        self.pwd_selected = Label(self.rootWindow,
                                  text = "None           "
                                  )
        self.pwd_selected.grid(row = 3,
                               column = 1,
                               columnspan = 3,
                               sticky = W
                               )


        output_label = Label(self.rootWindow,
                             text = "Ensemble output filename: "
                             )
        output_label.grid(row = 4,
                          column = 0,
                          sticky = W
                          )
        output = Entry(self.rootWindow,
                       textvariable = self.output
                       )
        output.grid(row = 4,
                    column = 1,
                    sticky=W
                    )


        gaps_label = Label(self.rootWindow,
                           text = "Chain-breaks permitted? "
                           )
        gaps_label.grid(row = 5,
                        column = 0,
                        sticky = W
                        )
        gaps_dropdown = OptionMenu(self.rootWindow,
                                   self.gap_setting,
                                   "None",
                                   "Some",
                                   "All",
                                   command=self.Semipermissive
                                   )
        gaps_dropdown.grid(row = 5,
                           column = 1,
                           sticky = E
                           )

        break_num_label = Label(self.rootWindow,
                                text = "# of chain-breaks permitted: "
                                )
        break_num_label.grid(row = 6,
                             column = 0,
                             sticky = W
                             )
        self.break_num_entry = Entry(self.rootWindow,
                                     textvariable = self.break_num
                                     )
        self.break_num_entry.grid(row = 6,
                                  column = 1,
                                  sticky=W
                                  )
        self.break_num_entry.configure(state='disable')


        align_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Perform sequence alignment",
                                        variable = self.align,
                                        command = self.align_check
                                        )
        align_checkbutton.grid(row = 7,
                               column = 0,
                               sticky = W+E
                               )




        self.template_entry = Button(self.rootWindow,
                                     text = "Select Template File",
                                     command = self.select_file
                                     )
        self.template_entry.grid(row = 8,
                                  column = 0,
                                  sticky=E
                                  )
        self.template_entry.configure(state='disable')

        self.template_selected = Label(self.rootWindow,
                                    text = "None           "
                                    )
        self.template_selected.grid(row = 8,
                                 column = 1,
                                 columnspan = 3,
                                 sticky = W
                                 )


        chain_label = Label(self.rootWindow,
                                text = "Chain ID for template: "
                                )
        chain_label.grid(row = 9,
                             column = 0,
                             sticky = W
                             )
        self.chain_entry = Entry(self.rootWindow,
                                     textvariable = self.chain
                                     )
        self.chain_entry.grid(row = 9,
                                  column = 1,
                                  sticky=W
                                  )
        self.chain_entry.configure(state='disable')




        model_label = Label(self.rootWindow,
                                text = "Model ID for template: "
                                )
        model_label.grid(row = 10,
                             column = 0,
                             sticky = W
                             )
        self.model_entry = Entry(self.rootWindow,
                                 text = "0",
                                 textvariable = self.model
                                     )
        self.model_entry.grid(row = 10,
                                  column = 1,
                                  sticky=W
                                  )
        self.model_entry.configure(state='disable')

        percent_label = Label(self.rootWindow,
                                text = "Percent ID Cutoff: "
                                )
        percent_label.grid(row = 11,
                             column = 0,
                             sticky = W
                             )
        self.percent_entry = Entry(self.rootWindow,
                                   text = "0.7",
                                   textvariable = self.percent
                                   )
        self.percent_entry.grid(row = 11,
                                  column = 1,
                                  sticky=W
                                  )
        self.percent_entry.configure(state='disable')



        go_button = Button(self.rootWindow, text = "Go!", command = self.execute);
        go_button.grid(row=12, column = 0, columnspan = 2, sticky = E+W);

        status_label = Label(self.rootWindow, text = "Status: ");
        status_label.grid(row = 13, column = 0, sticky = W);
        self.status = Label(self.rootWindow, text = "Waiting");
        self.status.grid(row = 13, column = 1, columnspan = 3, sticky = W);


    def align_check(self):
        if self.align.get() == 1:
            self.template_entry.configure(state = 'normal')
            self.chain_entry.configure(state = 'normal')
            self.model_entry.configure(state = 'normal')
            self.percent_entry.configure(state = 'normal')
        else:
            self.template_entry.configure(state = 'disable')
            self.chain_entry.configure(state = 'disable')
            self.model_entry.configure(state = 'disable')
            self.percent_entry.configure(state = 'disable')

    def Semipermissive(self, event):
        if self.gap_setting.get() == "Some":
            self.break_num_entry.configure(state='normal')
        else:
            self.break_num_entry.configure(state='disable')


    def select_files(self):
        self.status["text"] = "Waiting";
        self.inputfiles = self.utility.select_input_files(1,self.files_selected)

    def select_file(self):
        self.status["text"] = "Waiting";
        self.template = self.utility.select_input_file(self.template_selected)

    def select_pwd(self):
        self.status["text"] = "Waiting";
        self.dir_to_save = self.utility.select_output_dir(self.pwd_selected)

    def execute(self):
        self.status["text"] = "Waiting";

        break_num = int(self.break_num.get());
        gap_setting = self.gap_setting.get();



        ## Manually specify the arguments!

        # shared arguments
        options.pwd = self.dir_to_save

        # prepare input arguments
        options.input = self.inputfiles
        options.output = str(self.output.get())
        if gap_setting == "None":
            options.permissive = False
            options.semipermissive = 0
        elif gap_setting == "All":
            options.permissive = True
            options.semipermissive = 0
        elif gap_setting == "Some":
            options.permissive = False
            options.semipermissive = break_num
        else:
            options.permissive = False
            options.semipermissive = 0


        if self.align.get() == 1:
            options.align = True
            options.template = os.path.basename(self.template)\
                [0:(len(os.path.basename(self.template)) - 4)]

            options.chain = self.chain.get()
            options.model = self.model.get()
            options.percent = float(self.percent.get())

        else:
            options.align = False



        if(self.dir_to_save != "" and len(self.inputfiles) > 0):
            if (self.align.get() == 1 and self.template != "") \
                or (self.align.get() == 0):
                self.utility.prepare_input_command_run(options,
                                                       self.status,
                                                       self.rootWindow
                                                       )
                self.rootWindow.destroy()
            else:
                print("\n\nIf you have chosen to do an alignment, you must"
                      " also select a template file to use as your guide."
                      " Models that are less similar than the percent identity"
                      " cutoff will not be included in the ensemble."
                      )


        else:
            print("\n\nPlease select a working directory (to save results in)"
                  ", and some input files to make an ensemble out of!"
                  )




class Analyze:
    def __init__(self, utility):
        self.utility = utility;
        self.rootWindow = Tk();
        self.rootWindow.wm_title("Analyze Ensemble");
        self.dcut = StringVar(self.rootWindow);
        self.dcut.set("2.5");
        self.ensemble = ""
        self.groupm = StringVar(self.rootWindow)
        self.groupm.set("")
        self.groupn = StringVar(self.rootWindow)
        self.groupn.set("")
        self.maxclust = IntVar(self.rootWindow)
        self.maxclust.set(3)
        self.auto = IntVar(self.rootWindow)
        self.auto.set(0)
        self.avg = IntVar(self.rootWindow)
        self.avg.set(0)
        self.color = IntVar(self.rootWindow)
        self.color.set(0)
        self.dir_to_save = "";
        self.setup_gui();
        self.rootWindow.mainloop();

    def setup_gui(self):


        pwd_button = Button(self.rootWindow,
                            text = "Select Working Directory",
                            command = self.select_pwd
                            )
        pwd_button.grid(row=0,
                        column = 0,
                        columnspan = 2,
                        sticky = E+W
                        )

        pwd_selected_label = Label(self.rootWindow,
                                   text = "Working Directory: "
                                   )
        pwd_selected_label.grid(row = 1,
                                column = 0,
                                sticky = W
                                )
        self.pwd_selected = Label(self.rootWindow,
                                  text = "None           "
                                  )
        self.pwd_selected.grid(row = 1,
                               column = 1,
                               columnspan = 3,
                               sticky = W
                               )



        select_button = Button(self.rootWindow,
                               text = "Select Input Ensemble",
                               command = self.select_ensemble
                               )
        select_button.grid(row=2,
                           column = 0,
                           columnspan = 2,
                           sticky = E+W
                           );

        file_selected_label = Label(self.rootWindow,
                                     text = "Ensemble Selected: "
                                     )
        file_selected_label.grid(row = 3,
                                  column = 0,
                                  sticky = W
                                  )
        self.ensemble_selected = Label(self.rootWindow,
                                    text = "None           "
                                    )
        self.ensemble_selected.grid(row = 3,
                                 column = 1,
                                 columnspan = 3,
                                 sticky = W
                                 )



        dcut_label = Label(self.rootWindow,
                             text = "Cutoff distance for core atoms: "
                             )
        dcut_label.grid(row = 4,
                          column = 0,
                          sticky = W
                          )
        dcut = Entry(self.rootWindow,
                       textvariable = self.dcut
                       )
        dcut.grid(row = 4,
                    column = 1,
                    sticky=W
                    )


        groupm_label = Label(self.rootWindow,
                                text = "Group M models: "
                                )
        groupm_label.grid(row = 5,
                             column = 0,
                             sticky = W
                             )
        self.groupm_entry = Entry(self.rootWindow,
                                     textvariable = self.groupm
                                     )
        self.groupm_entry.grid(row = 5,
                                  column = 1,
                                  sticky=W
                                  )


        groupn_label = Label(self.rootWindow,
                                text = "Group N models: "
                                )
        groupn_label.grid(row = 6,
                             column = 0,
                             sticky = W
                             )
        self.groupn_entry = Entry(self.rootWindow,
                                 text = "0",
                                 textvariable = self.groupn
                                     )
        self.groupn_entry.grid(row = 6,
                                  column = 1,
                                  sticky=W
                                  )






        auto_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Perform automatic clustering",
                                        variable = self.auto,
                                        command = self.auto_check
                                        )
        auto_checkbutton.grid(row = 7,
                               column = 0,
                               columnspan = 2,
                               sticky = W
                               )



        maxclust_label = Label(self.rootWindow,
                                text = "Max # of clusters to search for: "
                                )
        maxclust_label.grid(row = 9,
                             column = 0,
                             sticky = W
                             )
        self.maxclust_entry = Entry(self.rootWindow,
                                     textvariable = self.maxclust
                                     )
        self.maxclust_entry.grid(row = 9,
                                  column = 1,
                                  sticky=W
                                  )
        self.maxclust_entry.configure(state='disable')




        avg_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Use average deviation rather than RMSD.",
                                        variable = self.avg
                                        )
        avg_checkbutton.grid(row = 10,
                               column = 0,
                               columnspan = 2,
                               sticky = W
                               )

        color_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Set b-factors in final " +\
                                               "ensemble equal to inter-LODR" +\
                                               " (or group M LODR).",
                                        variable = self.color
                                        )
        color_checkbutton.grid(row = 11,
                               column = 0,
                               columnspan = 2,
                               sticky = W
                               )




        go_button_analyze = Button(self.rootWindow,
                                   text = "Analyze!",
                                   command = self.execute_analyze
                                   )
        go_button_analyze.grid(row=12,
                               column = 0,
                               columnspan = 2,
                               sticky = E+W
                               )

        status_label = Label(self.rootWindow, text = "Status: ");
        status_label.grid(row = 13, column = 0, sticky = W);
        self.status = Label(self.rootWindow, text = "Waiting");
        self.status.grid(row = 13, column = 1, columnspan = 3, sticky = W);

    def auto_check(self):
        if self.auto.get() == 1:
            self.groupm_entry.configure(state = 'disable')
            self.groupn_entry.configure(state = 'disable')
            self.maxclust_entry.configure(state='normal')
        else:
            self.groupm_entry.configure(state = 'normal')
            self.groupn_entry.configure(state = 'normal')
            self.maxclust_entry.configure(state='disable')

    def select_ensemble(self):
        self.status["text"] = "Waiting";
        self.ensemble = self.utility.select_input_ensemble(self.ensemble_selected)

    def select_pwd(self):
        self.status["text"] = "Waiting";
        self.dir_to_save = self.utility.select_output_dir(self.pwd_selected)

    # parser for group numbers
    def parse_range(self, value):
        result = set()
        for part in value.split(','):
            x = part.split('-')
            result.update(range(int(x[0]), int(x[-1]) + 1))
        return(sorted(result))


    def execute_analyze(self):
        self.status["text"] = "Waiting";

        ## Manually specify the arguments!
        options.pwd = self.dir_to_save


        # analyze arguments
        options.input = self.ensemble
        options.dcut = float(self.dcut.get())


        options.auto = self.auto.get()
        if options.auto == 1:
            options.maxclust = self.maxclust.get()
        else:
            try:
                options.groupm = self.parse_range(str(self.groupm.get()))
            except:
                print "\n\nPlease select a group M, or choose the auto option."
                print "\nGroups should be formatted like so: 0,3,5,6-9,13\n\n"

            try:
                options.groupn = self.parse_range(str(self.groupn.get()))
            except:
                pass

        options.avg = True if (self.avg.get() == 1) else False
        options.color = True if (self.color.get() == 1) else False



        if(self.dir_to_save != "" and self.ensemble != "") and\
           (self.auto.get() == 1 or self.groupm.get() != ""):
                self.utility.analyze_command_run(options,
                                                       self.status,
                                                       self.rootWindow
                                                       )
                self.rootWindow.destroy()
        else:
            print("\n\nPlease select a working directory (to save results in)"
                  ", a prepared ensemble to analyze, and either select the "
                  "auto-cluster option, or define at least a group M to "
                  "analyze. Then click the button again!"
                  )




main_window = MainRoot();
