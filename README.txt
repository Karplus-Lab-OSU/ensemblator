NOTE: This readme is out of date as of now (Feb. 28th 2016). 
      Since adding a GUI, the ENSEMBLATOR is now version 3. 
      I'm going to fix this readme soon.

ENSEMBLATOR version 2
Andrew E. Brereton - Oregon State University - 2015

1. Introduction:
	Hello! Welcome to the Ensemblator. The purpose of this readme is to help 
	you install and use the Ensemblator to best effect. The purpose of this 
	program is to create ensembles from pdb files, and then to intelligently 
	compare those ensembles. Many useful insights can be gained from directly 
	comparing two whole ensembles; as well, there is an option to automatically
	try to identify different groups within a single ensemble, and compare 
	those.
	
	This readme contains these parts:
	    1. Self-aware introduction!
	    2. Installation
	    3. Requirements
	    4. Usage
	        4a. prepare_input.py
	        4b. ensemblator.py
	    5. Warnings
	    6. Bug Reports
	    7. FAQ
	
2. Installation:
    To install the Ensemblator, simply ensure that "ensemblator.py" and 
    "prepare_input.py" are in your path, and are executable. Alternatively, you
    could place the scripts in the directory you wish to run them in, though
    that is much more tedious as you will have to copy them every time into a 
    new directory, or specify a long path.

3. Requirements:
    To use the Ensemblator you will require Python 2 to be installed on your
    system. Furthermore, you will need the following packages:
    
        numpy
        biopython
        matplotlib
    
    If you want to use the automatic clustering method to identify groups
    within your ensemble, you will also need to install the "scipy" package.
    This is recommended.
    
    In order to use the --align method of 'prepare_input.py', to align
    homologous or similar structures by sequence, you will need the alignment
    program 'Muscle' to be installed on your computer.

    
4. Usage:
    Both "ensemblator.py" and "prepare_input.py" have a list of all 
    command-line options that is accessible using '-h'. Before an ensemble can
    be analyzed with the Ensemblator, the input must be prepared by using the
    "prepare_input.py" script.
    
    4a. prepare_input.py:
        The purpose of this script is to compare a series of structures, and
        create an ensemble which only contains all the atoms common to each of
        the structures. It does this using the residue ID and the atom type.
        ie. All the residue 46 CA atoms need to be present in all structures to
        be included. This is important to consider, as if the numbering of the
        residues is different, only backbone atoms will be kept, and the
        subsequent analysis by the Ensemblator will be weird.
        There are always at least two outputs from this script: your output 
        ensemble, and a legend that will contain the model id associated with
        each input structure. 
    *** It is important to view this legend, as the nature of the algorithm
    *** for assembling the ensemble means that the models will be in a random
    *** order in the final ensemble.
        
        -h
            Shows a help dialogue.
        -i or --input
            Used to specify the input. This should be a pdb file with at least
            two models, chains, or alternate conformations. Furthermore, any
            number of pdb files can be used, and all will be combined into
            the final prepared ensemble.
        -o or --output
            Used to specify the final name of the output ensemble.
        -p or --permissive
            By default, the "prepare_input.py" script is not permissive. What
            this means is that if the script detects that there is a 
            chain-break or a gap in the backbone of the structure, it will not
            include this structure in the final analysis. However, if you run
            the script with this option, it will include these structures in
            the final analysis, meaning that all structures will have a gap at
            that position.
        --semipermissive
            As above in the case of permissive, except here, it will still
            remove structures if they are missing more than a specified number
            of residues. ie. "--semipermissive 3" will remove any structure 
            with more than three gaps in any positions. This option will also
            remove structures that have many ambiguous missing atoms. At one
            point, the program lists all backbone atoms, then breaks them
            into sets of four, in the order they appear in the file. If there
            is only one unique set (N,CA,C,O), it means that there are no
            residues that are missing backbone atoms. If there are 2 sets, it
            means there is at least one residue that is missing backbone atoms.
            With this check in place, it will make sure that the number sets
            does not exceed the number given by this option + 1. ie. If
            "--semipermissive 1" it will keep structures with only one gap, and
            having two or less types of backbone order. HOWEVER! For this atom
            order check, there is a maximum number of types of backbone order
            of 6. ie, no matter what value of semipermissive is set, it will
            still remove structures with more than 6 types of backbone atom
            order. Anything more than that is a seriously messed up pdb file,
            and of questionable use and trustworthiness.
        -s or --skipxray
            This will skip the initial formatting step that runs on the pdb
            files, which separates them into model, chain, and altconf files.
            This should only ever be used if you have a series of identical
            structures, with only ATOM lines, with only one model, chain, and
            alt conf per file. This will save some computation time, but is
            only really noticeable when dealing with a very very large number
            of files.
        -l or --log
            Setting this option will have the script save a log to the working
            directory with more details about what it was doing.
        -a or --align
            Setting this option will do multiple sequence alignment using 
            MUSCLE, which must be installed on the user's computer. This
            multiple sequence alignment will be used to renumber the structures
            to ensure that all residues are correctly matched.
            This will output a file: 'muscle_align.fasta'. This file is the
            fasta formatted alignment that was used to generate the new
            numbering of the residues.
        -c or --chain
            This option is required when using --align. This option is used to
            define a template that all the aligned sequences will be compared
            to. If they are less than a certain percent identity (by default
            70%) they will not be included in the analysis. The format for this
            argument is very important. 
            Here are some examples:
                --chain 3fc2.pdb,A
                    This will select chain A, model 0, from the file 3fc2.pdb
                -c 3fc3.pdb,B,20
                    This will select chain B, model 20, from the file 3fc2.pdb
                --chain 3fc3.pdb
                    THIS WILL NOT WORK. At least a chain and file must be
                    specified.
                -c 3fc2.pdb,0
                    THIS WILL NOT WORK. This will try to use "0" as the chain
                    ID, which won't work. The chain ID must ALWAYS be 
                    specified.
        --percent
            This will allow you to set the percent identity to use when
            deciding to keep structures or not. You can use percents (ie. 70, 
            100, or 35.6546), or you can use fractional percents (ie. 0.5 will
            be treated as 50%).
    
    4b. ensemblator.py:
        This program is what does the actual analysis of your ensemble. It will
        output a few tables and graphs, depending on the settings you run it
        with:
            'pairwise_analysis.tsv':
                This tab-separated table contains information about each of the
                pairs of models. From left to right, the columns list: the id 
                of model x, the id of model y, the number of atoms removed from
                the core for this pair, the rmsd for all the atoms in the
                two structures, and the rmsd for only the core atoms in this
                pair of structures.
            'eeGlobal_out.tsv':
                This tab-separated table contains information each atom in the
                ensemble. From left to right the columns describe: the residue
                id of the atom, the atom type, the RMSD of the atom calculated
                pairwise from group M (ie. the RMSD of all the pairwise
                distances in group M), the same for group N, the same but
                calculated from each M to N pair, the closest distance between
                any member of M with any member of N, the pair of models which
                actually had that closest approach, and whether or not this 
                atom was included in the common core calculated for the 
                ensemble.
            'eeLocal_out.tsv':
                This tab-separated table contains information about the LODR
                calculated for each residue. The columns list from left to
                right: the residue id, the RMS of the LODR calculated for each
                pair of structures in group M, the same for group N, the same
                for each M to N pair, the minimum LODR for any member of M
                compared with any member of N, and which pair was that closest
                LODR score.
            'global_overlay_X.X.pdb':
                This overlay of structures is the overlay calculated by using
                the first model in the ensemble as a reference structure, and
                aligning all the other models to this first model, using only
                the common core atoms determined depending on your dcut score.
            'groups_legend.tsv':
                Only output if the --auto option is used.
                This tab-separated table contains information about which
                models are clustered into which groups. The left column
                contains the model id, and the right column the group that
                model was assigned to.
            'eeGLOBAL_dcut=X.X.png':
                A graph of some of the data from 'eeGlobal_out.tsv'.
            'eeLOCAL_dcut=X.X.png':
                A graph of some of the data from 'eeLocal_out.tsv'.
            'groupX.pdb':
                A pdb file containing only the members of group X. Will only
                be output if the --auto option is set.
        
        The options available to the user of this script are:
            -h
                Shows a help dialogue.
            -i or --input
                A pdb file that has been prepared by 'prepare_input.py'.
            -d or --dcut
                A value in Angstroms to use as a distance cutoff to define the
                common core of your ensemble. This can be as oddly specific as
                you'd like, so feel free to use 3.445352362362 if you so
                choose. The default value for this calculation is 2.5. Choosing
                this value is important, and different dcut values will often
                give different results. It is valuable to play around to
                determine what works best for your ensemble.
            --auto
                This option will allow the user to avoid telling the 
                Ensemblator which groups to compare. Instead, the program will
                do all the pairwise analysis, and then use these results to
                determine which statistics (# of atoms removed, rmsd of all 
                atoms, rmsd of core atoms) give the best clusters. There is a
                penalty for increasing numbers of clusters, which biases the
                discovery of clusters to lower numbers of clusters. Clustering
                is done using a k-means algorithm. The clustering algorithms
                will also disfavor a solution that has a cluster with only one
                member.
            --maxclust
                Allows the user to specify a maximum number of clusters to 
                identify within the ensemble. By default this number is 6. This
                can be increased as high as the user wants, or as low as 2.
                Higher maxclust numbers will slightly increase the computation
                time.
            -m or --groupm
                Define group M for analysis. If not using the auto option, then
                at least group M must be defined. Members of a group can be
                separated by commas, as well as ranges specified using dashes.
                For example, to specify all 20 members of an ensemble as group
                M, you would type '-m 0-19'. To specify only some, you might
                type '-m 0-4,13-19'.
            -n or --groupn
                Define group N for analysis and comparison to group M.
            -l or --log
                Setting this option will have the script save a log to the 
                working directory with more details about what it was doing.
            --skiplocal
                Setting this option will skip eeLocal analysis, ie, it will
                skip calculation of the LODR values. This can save time when
                changing only dcut, as if the groups are the same, the LODR
                scores will be identical to previous analysis, regardless of
                a new dcut value.
            --color
                Setting this will result in the final models output having the
                b factors replaced with the Inter-group (if more than one
                group) or Group M LODR. This allows easy visualization in 
                pymol using the "spectrum b" command.
                
5. Warnings:
    1.	All files used as input must be IN the present working directory. 
	    Otherwise, intermediate files are created in distant directories, and 
	    then cannot be access by the software.
	2.	If you want to keep all the files of a previous run, but wish to change
	    a variable and run the Ensemblator again, you will need to move or 
	    rename many of the output files. These files tend to always have the 
	    same name, and will be overwritten.
	3.	I have noticed that sometimes the "prepare_input.py" script does not 
	    work on pdb files directly generated by phenix. This has to do with the
	    names of the heavy atoms in the pdb files. If any of these are
	    lowercase, it will break the ensemblator.

6. Bug Reports:
	Please submit any issues if you have a bug!
	
