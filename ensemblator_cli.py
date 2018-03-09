#!/usr/bin/env python2


from argparse import ArgumentParser
from ensemblator.ensemblator_core import *


# checks if there is a valid file at a specified location
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)


parser = ArgumentParser()
parser.add_argument("--prepare",
                  action="store_true",
                  dest="prepare",
                  default=False,
                  help=("Use the ensemblator to create an ensemble ready"
                        " for analysis. Must be done at least once before"
                        " the analysis option can be used without erros."))
parser.add_argument("--analyze",
                  action="store_true",
                  dest="analyze",
                  default=False,
                  help=("Use the ensemblator to analyze an ensemble prepared"
                        " for analysis."))
parser.add_argument('-i',
                    '--input',
                    dest="input",
                    metavar=str,
                    help='This should be a pdb or cif file, or series of pdb or cif files.',
                    nargs='*')
parser.add_argument("-o",
                    "--output",
                    dest="output",
                    type=str,
                    help=("This descriptor will define the final name of the " + "output file."))
parser.add_argument("--pwd",
                  dest="pwd",
                  type=str,
                  help=("This defines a working directory to save all output files in."))
parser.add_argument("-p",
                  "--permissive",
                  action="store_true",
                  dest="permissive",
                  default=False,
                  help="If set, will use files to generate ensemble even" +
                       " if there are gaps or chainbreaks. These gaps will " +
                       "propigate into all members of the ensemble.")
parser.add_argument("--semipermissive",
                  dest="semipermissive",
                  type=int,
                  default=0,
                  help="If set, will use files to generate ensemble even" +
                       " if there are gaps or chainbreaks, but only if they " +
                       "contain less than the specified number of gaps"
                       " (missing residues). Will also remove structures if they"
                       " are too disordered (ie. have many ambigious missing "
                       " atoms.")
parser.add_argument(
    "-l",
    "--log",
    action="store_true",
    dest="log",
    default=False,
    help="If set, will generate a log file:" + " 'prepare_input.log'.")
parser.add_argument("--align",
                  action="store_true",
                  dest="align",
                  default=False,
                  help="If set, will perform a sequence alignment first,"
                       " and use those residue numbers to create a better result.")
parser.add_argument("-c",
                  "--chain",
                  type=str,
                  dest="template",
                  help="Will determine which structure to use when aligning"
                       " sequences. Select a specific input file, and then use a"
                       " comma and state which chain of the file to use as a"
                       " template for the alignments. See the README for"
                       " examples.")
parser.add_argument("--percent",
                  dest="percent",
                  default=0.7,
                  type=float,
                  help="Percent identity to use when building ensemble"
                       " using the -align function. Any structure with less"
                       " than this percent identity to the template chain"
                       " won't be included in the analysis.")
parser.add_argument("-d",
                  "--dcut",
                  dest="dcut",
                  type=float,
                  default=2.5,
                  help="Distance cutoff for core decision.")
parser.add_argument("--auto_cutoff",
                  dest="dcutAuto",
                  action="store_true",
                  default=False,
                  help="If set, will perform try to automatically detect a"
                       " distance cutoff that will result in a core with"
                       " between 20% and 40% of the total atoms.")
parser.add_argument(
    "--auto",
    action="store_true",
    dest="auto",
    default=False,
    help=("If set, will cluster the pairwise results and rerun "
          "the analysis using the detected groups. Will ignore "
          "manually set m and n values while used. In the final "
          "outputs, the first group in a filename corresponds to "
          "group M, and the second to N. "
          "ie. eeLocal_Group_0_Group_1.svg has group 0 as M and "
          "group 1 as N."))
parser.add_argument(
    "--maxclust",
    type=int,
    dest="maxclust",
    default=3,
    help=("Maximum number of clusters to group the results into."))
parser.add_argument(
    "--cores",
    type=int,
    dest="cores",
    default=4,
    help=("Number of cores to use for pairwise analysis"))
parser.add_argument("-g",
                    "--groups",
                    dest="groups",
                    type=str,
                    help=("Groups of which models to compare. Use "
                        "dashes for a range, commas separate entries. Spaces to separate multiple groups. "
                        "E.g. 1,3,5,8-19,23 2,4,6-7,20-22 24-30"),
                    nargs='*')
options = parser.parse_args()
# required options
if not options.prepare and not options.analyze:  # if filename is not given
    parser.error('Must choose either the --prepare or the --analyze option.')
if not options.pwd:
    parser.error('Must use the "--pwd" option to specify a working directory.')
if options.auto == False:
    if not options.groups and not options.prepare:
        parser.error('At least Group M must be specified, or turn on --auto')
if not options.input:  # if filename is not given
    parser.error('Input not specified')
if options.maxclust < 2:
    parser.error('Minimum of 2 clusters. Try again with a higher maxclust')
if not options.output and not options.analyze:  # if filename is not given
    parser.error('output filename not given')
if options.align:
    if not options.template:  # if filename is not given
        parser.error('Template not specified. When aligning, use '
                     '"--chain filename.pdb,X,X" to specifiy a template chain'
                     ' and model. eg. "--chain 1q4k.pdb,0,A"')
if options.analyze and len(options.input) > 1:
    parser.error('Only one input file may be specified for analysis!')
if not os.path.exists(options.pwd):
    os.makedirs(options.pwd)


if options.prepare == True and options.analyze == True:

    sys.exit("Please choose only one of the options"
             " '--prepare' or '--analyze'.")

elif options.prepare == False and options.analyze == False:

    sys.exit("Please choose only one of the options"
             " '--prepare' or '--analyze'.")


################################################################################
#### Prepare Input Files for Analysis ####
################################################################################

elif options.prepare == True and options.analyze == False:

    if options.align:

        template = []
        try:
            template.extend((options.template).split(","))
        except:
            sys.exit("Template must be formated like so:"
                     " '-chain filename.pdb,chain,model' e.g. "
                     " '-chain 1q4k.pdb,2,A' or '-chain 1q4k.pdb,A'.")

        try:
            try:
                options.template = template[0][0:len(template[0]) - 4]
                options.model = template[1]
                options.chain = template[2]

            # otherwise just use zero
            except:
                options.template = template[0][0:len(template[0]) - 4]
                options.model = 0
                options.chain = template[1]
        except:
            sys.exit("Template must be formated like so:"
                     " '-chain filename.pdb,chain,model' e.g. "
                     " '-chain 1q4k.pdb,2,A' or '-chain 1q4k.pdb,A'.")

    prepare_input(options)


################################################################################
#### Analyze ####
################################################################################

elif options.analyze == True and options.prepare == False:

    analyze(options)
