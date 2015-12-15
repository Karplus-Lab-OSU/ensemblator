#!/usr/bin/env python

import numpy
import time
import fnmatch
import Bio
from optparse import OptionParser
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import os
import sys
import re
from itertools import permutations

# function for iterating
def grouped(iterable, n):
    # s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), 
    # (s2n,s2n+1,s2n+2,...s3n-1), ...
    return izip(*[iter(iterable)]*n)

#function to allow multiple pdb files in the command line
def cb(option, opt_str, value, parser):
    args=[]
    for arg in parser.rargs:
        if arg[0] != "-":
            args.append(arg)
        else:
            del parser.rargs[:len(args)]
            break
        if getattr(parser.values, option.dest):
            args.extend(getattr(parser.values, option.dest))
    setattr(parser.values, option.dest, args)
# parser for command line options
parser=OptionParser()
parser.add_option(
                    "-i", 
                    "--input", 
                    dest="input",
                    action="callback", 
                    callback=cb, 
                    metavar="FILE", 
                    help="This should be a pdb file, or series of pdb files."
                    )
parser.add_option(
                    "-o", 
                    "--output", 
                    dest="output",
                    type="str", 
                    help="This descriptor will define the final name of the " +
                    "output file."
                    )
parser.add_option(
                    "-p", 
                    "--permissive", 
                    action="store_true", 
                    dest="permissive", 
                    default=False, 
                    help="If set, will use files to generate ensemble even" +
                    " if there are gaps or chainbreaks. These gaps will " +
                    "propigate into all members of the ensemble."
                    )
parser.add_option(
                    "--semipermissive",  
                    dest="semipermissive",
                    type="int", 
                    default=0, 
                    help="If set, will use files to generate ensemble even" +
                    " if there are gaps or chainbreaks, but only if they " +
                    "contain less than the specified number of gaps"
                    " (missing residues). Will also remove structures if they"
                    " are too disordered (ie. have many ambigious missing "
                    " atoms."
                    )
parser.add_option(
                    "-s", "--skipxray", 
                    action="store_true", 
                    dest="skip", 
                    default=False, 
                    help="If set, will skip xray-prep module. Only use if" + 
                    " you know each file contains one model, one chain, " + 
                    "and one conformation."
                    )
parser.add_option(
                    "-l", 
                    "--log", 
                    action="store_true", 
                    dest="log", 
                    default=False, 
                    help="If set, will generate a log file:" + 
                    " 'prepare_input.log'."
                    )
parser.add_option(
                    "-a", 
                    "--align", 
                    action="store_true", 
                    dest="align", 
                    default=False, 
                    help="If set, will performa a sequence alignment first,"
                    " and use those residue numbers to create a better result."
                    )
parser.add_option(
                    "-c", 
                    "--chain", 
                    type="str",
                    dest="template", 
                    help="Will determine which structure to use when aligning"
                    " sequences. Select a specific input file, and then use a"
                    " comma and state which chain of the file to use as a"
                    " template for the alignments. See the README for"
                    " examples."
                    )
parser.add_option(
                    "--percent",  
                    dest="percent", 
                    default=0.7,
                    type="float",
                    help="Percent identity to use when building ensemble"
                    " using the -align function. Any structure with less"
                    " than this percent identity to the template chain"
                    " won't be included in the analysis."
                    )
(options, args) = parser.parse_args()

if not options.output:   # if filename is not given
    parser.error('output filename not given')
if not options.input:   # if filename is not given
    parser.error('Input not specified')
if options.align:
    if not options.template:   # if filename is not given
        parser.error('Template not specified. When aligning, use '
                     '"--chain filename.pdb,X,X" to specifiy a template chain'
                     ' and model. eg. "--chain 1q4k.pdb,A,0"')


# writes a log file
log = open("prepare_input.log", "a")



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
                numpy.abs(super_imposer.rms) < 0.0000001
                assert \
                numpy.max(numpy.abs(super_imposer.rotran[1])) < 0.000001
                assert \
                numpy.max(numpy.abs(super_imposer.rotran[0]) - \
                                            numpy.identity(3)) < 0.000001
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
    os.remove(pdb)
    fileout = open(pdb, 'w')
    
    for line in filein:
        if line[0:3] == 'END' and line[0:6] != 'ENDMDL':
            pass
        else:
            fileout.write(line)
    fileout.write("END")
    fileout.close()
    filein.close()


# eeprep function, is passed a list of xray-prepped files
# It's purpose is to check for equvilence of each atom type at each 
# residue in all the structures
# it removes non equivelent atoms from the ensemble
def eeprep(pdbs):

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
    n = max(resnum_count_list)
        
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
    for key in structures:
        io.set_structure(structures[key])
        io.save(str(key)+'_out.pdb')

    outputname = options.output

    outfile = open(outputname, 'w')
    outfile.write("")
    outfile.close

    # formatting the output files here
    counter = 0
    for key in structures:
        filename = str(key)+'_out.pdb'                        
        if backbone_scan(filename) == True:
            os.remove(filename)
            legend_dict[key] = 'REMOVED FROM FINAL ENSEMBLE'
        else:
            infile = open(filename,'r')
            os.remove(filename)
            outfile = open(outputname, 'a')
            # formatting to ensure that the MODEL line is correct
            modeltag = str(counter)
            while len(modeltag) < 9:
                modeltag = " " + modeltag
            if len(modeltag) == 9:
                modeltag = "MODEL" + modeltag + "\n" 
            # write model line
            outfile.write(modeltag)
            # write all the atom lines
            for line in infile:
                if line[0:6] == 'ATOM  ':
                    outfile.write(line)
            # write endmdl line
            outfile.write("ENDMDL\n")
            outfile.close()
            counter += 1
    # need to cap the file with 'END'
    outfile = open(outputname, 'a')
    outfile.write("END   \n")
    outfile.close()

    # rewrites the files to correctly format them, removing the alt_conf id    
    infile = open(outputname, 'r')
    os.remove(outputname)
    outfile = open(outputname, 'w')               
    for line in infile:    
        if line[0:6] == 'ATOM  ':
            line = line[0:16] + " " + line[17:len(line)]
            outfile.write(line)
        else:
            outfile.write(line)
    return legend_dict

# xray-prep function, reads in a single pdb file at a time, 
# also needs an output name to iterate over
# many variable names in Esperanto, sorry
def xrayprep(pdb, output):

    pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)

    # reads the pdb file
    try:
        strukturo = pdb_reader.get_structure("temp", pdb)
    except:
        sys.exit("\n\n\nFailed to read file:"
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
        endosiero = open(outputname,'r')
        os.remove(outputname)
        eldosiero = open(outputname, 'w')
    

        for line in endosiero:
            if line[0:6] == 'ATOM  ':
                new_line = line[0:16] +\
                             " " + \
                             line[17:21] + \
                             "A" + \
                             line[22:len(line)]
                eldosiero.write(new_line)
            elif line[0:6] == 'HETATM':
                pass
            else:
                new_line = line
                eldosiero.write(new_line)
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
          
    endosiero = open(outputname,'r')
    os.remove(outputname)
    eldosiero = open(outputname, 'w')
                    
    eldosiero.write("MODEL        0\n")
                
                
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
                eldosiero.write(new_line + otherStore)
                nPrint = False
                aPrint = False
                cPrint = False
                oPrint = False   
                                               
                prev_resnum = resnum
                otherStore = ""   
                new_line = "" 
                                    
                                    
    eldosiero.write("TER   \n")        
    eldosiero.write("ENDMDL")

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
        log.write("No ATOM lines in model " + 
                  str(pdb)[0:(len(str(pdb))-8)] + 
                  ", removing from final ensemble.\n"
                  )
        return True





#############################################################################
#############################################################################
#############################################################################
# mark start time

startTime = time.time()        
# Here the real code begins:

pdbs = options.input                                        
prepped_files_list = []


# skips xrayprep
if options.skip == True:
    for pdb in pdbs:
        print "Formatting " + str(pdb) + " for input."
        log.write( "Formatting " + str(pdb) + " for input.\n")    
        xray_outname = "prepped_" + str(pdb[0:(len(pdb)-4)])
        filein = open(pdb, 'r')
        fileout = open(xray_outname, 'w')
        fileout.write("MODEL        0\n")
        for line in filein:
            if line[0:6] == 'ATOM  ':
                fileout.write(line)
        fileout.write("TER   \n")        
        fileout.write("ENDMDL")
        fileout.close()
# runs xray-prep
else:
    print("\n##################\n\nFormatting and seperating..."
          "\n\n##################\n"
          )
    log.write("\n##################\n\nFormatting and seperating..."
              "\n\n##################\n"
              )
    for pdb in pdbs:
        print("Formatting " + 
              str(pdb) + 
              " for input, and seperating into " +
              "model, chain, and altconf files."
              )
        log.write("Formatting " + 
                  str(pdb) + 
                  " for input, and seperating into" +
                  " model, chain, and altconf files.\n"
                  )
        xray_outname = "prepped_" + str(pdb[0:(len(pdb)-4)])
        prepped_files_list.extend(xrayprep(pdb, xray_outname))



# runs backbone_scan
if options.permissive == False:
    print("\n##################\n\nScanning backbones for gaps and "
          "missing residues...\nTo prevent this, use option"
          " '--permissive'.\n\n##################\n"
          )
    log.write("\n##################\n\nScanning backbones for gaps"
              " and missing residues...\nTo prevent this, use option"
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
        log.write ("Removing file: " + 
                   str(prepped_file) + 
                   " from analysis, due to a gap,\nor atoms being in the" +
                   " wrong order in the file (this breaks the " +
                   "ensemblator).\n\n"
                   )
        bad_files.append(prepped_file)




good_files = set(prepped_files_list) - set(bad_files)
good_files = list(good_files)


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
                    pass
            seqs.write(">" + str(pdb) + "\n" + str(seq) + "\n")
        seqs.close()

        # command line args for MUSCLE
        cline = MuscleCommandline(input="sequences.fasta",
                                  out="muscle_align.fasta"
                                  )    
        log_str = "\n##################\n\n" + \
                  "Running MUSCLE alignment program with the following " +\
                  "command: \n" + str(cline) +\
                  "\n\n##################\n"
        print log_str
        log.write(log_str)
        
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
        

        template = []
        try:
            template.extend((options.template).split(","))
        except:
            sys.exit("Template must be formated like so:"
                     " '-chain filename.pdb,chain,model' e.g. "
                     " '-chain 1q4k.pdb,A,2' or '-chain 1q4k.pdb,A'."
                    )

        
        try:
            try:
                template_name = "prepped_" + \
                                template[0]\
                                    [0:(len(template[0]) - 4)] + \
                                "_model_" + \
                                str(template[2]) + \
                                "_chain_" + \
                                str(template[1]) + \
                                "_alt_A.pdb"
            # otherwise just use zero
            except:    
                template_name = "prepped_" + \
                                template[0]\
                                    [0:(len(template[0]) - 4)] + \
                                "_model_0" + \
                                "_chain_" + \
                                str(template[1]) + \
                                "_alt_A.pdb"
        except:
            sys.exit("Template must be formated like so:"
                     " '-chain filename.pdb,chain,model' e.g. "
                     " '-chain 1q4k.pdb,A,2' or '-chain 1q4k.pdb,A'."
                    ) 
        
        
        
        #empty the good_files list, to refill with the newly good files
        # ie the files with a certain percent id
        good_files = []
        align_counter = 0                 
        for key in aligned_dict:
            try:
                template = aligned_dict[template_name]
            except:
                sys.exit("Couldn't load template file, probably it was removed"
                         " after backbone scanning. Try running with the"
                         " --permissive option on, or pick a different "
                         "template."
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
                log.write(message + "\n")
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
log.write("\n##################\n\nCombining all files into an ensemble " +
          "ready for use with ensemblator by ensuring that all atoms match" +
          " in all structures...\n\n##################\n"
          )

    

# runs eeprep
legend_dict = eeprep(good_files)
legend_file = open('prepared_input_legend.tsv', 'w')
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
log.write("\n##################\n\nWrote file: " + 
          str(options.output) + 
          "\n\n##################\n"
          )


# align all the structures to make the output look nicer

print "\n##################\n\nCalculating overlay...\n\n##################\n"
log.write("\n##################\n\n"
          "Calculating overlay..."
          "\n\n##################\n"
          )
try:
    aligner(options.output)
except:
    err_str=("The overlay failed. This usually means that there were no"
          " atoms in the ensemble, which in turn usually means that afer"
          " removing the atoms not common to all structures, you were left"
          " with nothing. Sometimes, a single chain will simply be a few"
          " continious residues with no gaps, this can be greatly disruptive."
          " Try rerunning with --permissive set."
          " Or, even better, remove these bad chains from the analysis.\n\n\n"
          )
    print(err_str)
    log.write(err_str)

    
    
# time stuff
endTime = time.time()
workTime =  endTime - startTime
print "\n\nCompleted in: " + str(workTime) + " seconds.\n"
log.write( "\n\nCompleted in: " + str(workTime) + " seconds.\n")

log.close()

# removes the log file if the user didn't ask for it
if options.log == True:
    pass
else:
    os.remove("prepare_input.log")
