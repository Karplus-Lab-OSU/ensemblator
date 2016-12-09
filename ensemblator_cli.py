#!/usr/bin/env python2


from optparse import OptionParser
from ensemblator_core import *


# function to allow multiple pdb files in the command line
def cb(option, opt_str, value, parser):
    args = []
    for arg in parser.rargs:
        if arg[0] != "-":
            args.append(arg)
        else:
            del parser.rargs[:len(args)]
            break
        if getattr(parser.values, option.dest):
            args.extend(getattr(parser.values, option.dest))
    setattr(parser.values, option.dest, args)


parser = OptionParser()
parser.add_option("--prepare",
                  action="store_true",
                  dest="prepare",
                  default=False,
                  help=("Use the ensemblator to create an ensemble ready"
                        " for analysis. Must be done at least once before"
                        " the analysis option can be used without erros."))
parser.add_option("--analyze",
                  action="store_true",
                  dest="analyze",
                  default=False,
                  help=("Use the ensemblator to analyze an ensemble prepared"
                        " for analysis."))
parser.add_option("-i",
                  "--input",
                  dest="input",
                  action="callback",
                  callback=cb,
                  metavar="FILE",
                  help="This should be a pdb file, or series of pdb files.")
parser.add_option(
    "-o",
    "--output",
    dest="output",
    type="str",
    help="This descriptor will define the final name of the " + "output file.")
parser.add_option("--pwd",
                  dest="pwd",
                  type="str",
                  help="This defines a working directory to save all output"
                       " files in.")
parser.add_option("-p",
                  "--permissive",
                  action="store_true",
                  dest="permissive",
                  default=False,
                  help="If set, will use files to generate ensemble even" +
                       " if there are gaps or chainbreaks. These gaps will " +
                       "propigate into all members of the ensemble.")
parser.add_option("--semipermissive",
                  dest="semipermissive",
                  type="int",
                  default=0,
                  help="If set, will use files to generate ensemble even" +
                       " if there are gaps or chainbreaks, but only if they " +
                       "contain less than the specified number of gaps"
                       " (missing residues). Will also remove structures if they"
                       " are too disordered (ie. have many ambigious missing "
                       " atoms.")
parser.add_option(
    "-l",
    "--log",
    action="store_true",
    dest="log",
    default=False,
    help="If set, will generate a log file:" + " 'prepare_input.log'.")
parser.add_option("--align",
                  action="store_true",
                  dest="align",
                  default=False,
                  help="If set, will performa a sequence alignment first,"
                       " and use those residue numbers to create a better result.")
parser.add_option("-c",
                  "--chain",
                  type="str",
                  dest="template",
                  help="Will determine which structure to use when aligning"
                       " sequences. Select a specific input file, and then use a"
                       " comma and state which chain of the file to use as a"
                       " template for the alignments. See the README for"
                       " examples.")
parser.add_option("--percent",
                  dest="percent",
                  default=0.7,
                  type="float",
                  help="Percent identity to use when building ensemble"
                       " using the -align function. Any structure with less"
                       " than this percent identity to the template chain"
                       " won't be included in the analysis.")
parser.add_option("-d",
                  "--dcut",
                  dest="dcut",
                  type='float',
                  default=2.5,
                  help="Distance cutoff for core decision.")
parser.add_option(
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
parser.add_option(
    "--maxclust",
    type='int',
    dest="maxclust",
    default=3,
    help=("Maximum number of clusters to group the results into."))
parser.add_option("-m",
                  "--groupm",
                  dest="groupm",
                  type='str',
                  action="callback",
                  callback=parse_range,
                  help=("Models to use in comparision in the first group. Use "
                        "dashes for a range, commas separate entries. "
                        "E.g. 1,3,5,8-19,23"))
parser.add_option("-n",
                  "--groupn",
                  dest="groupn",
                  type='str',
                  action="callback",
                  callback=parse_range,
                  help=("Models to use in comparision in the second group. "
                        "Use dashes for a range, commas separate entries. "
                        "E.g. 1,3,5,8-19,23. OPTIONAL: do not include a "
                        "group N in order to do a single ensemble analysis."))
parser.add_option("--color",
                  action="store_true",
                  dest="color",
                  default=False,
                  help=("If set, will set b-factors in the final overlay to "
                        "relative inter-group LODR scores (if possible, "
                        "otherwise uses group m scores). Will not do "
                        "anything in auto-mode."))
(options, args) = parser.parse_args()
# required options
if not options.prepare and not options.analyze:  # if filename is not given
    parser.error('Must choose either the --prepare or the --analyze option.')
if not options.pwd:
    parser.error('Must use the "--pwd" option to specify a working directory.')
if options.auto == False:
    if not options.groupm and not options.prepare:
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
                     ' and model. eg. "--chain 1q4k.pdb,A,0"')
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
                     " '-chain 1q4k.pdb,A,2' or '-chain 1q4k.pdb,A'.")

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
                     " '-chain 1q4k.pdb,A,2' or '-chain 1q4k.pdb,A'.")

    prepare_input(options)


################################################################################
#### Analyze ####
################################################################################

elif options.analyze == True and options.prepare == False:

    analyze(options)
