ENSEMBLATOR version 3
=====================

Andrew E. Brereton - Oregon State University - 2016
---------------------------------------------------

Hello! Welcome to the *Ensemblator*. The purpose of this readme is to help you install and use the *Ensemblator* to best effect. The purpose of this program is to create ensembles from pdb files, and then to intelligently compare those ensembles. Many useful insights can be gained from directly comparing two whole ensembles; as well, there is an option to automatically try to identify different groups within a single ensemble, and compare those.

![image]

### Table of Contents

-   [Installation]
-   [Requirements]
-   [Usage]

    > -   [The main menu]
    > -   [Preparing your files for analysis]
    > -   [Analyzing your prepared ensemble]
    > -   [Understanding the output]

-   [Known Bugs]
-   [Bug Reports]

### Installation:

There are a few ways to install and use the *Ensemblator*:

1.  Run the source code:
    1.  [ensemblator.py] - The GUI version of the *Ensemblator*. As above you will need to make it executable, or run it using python from the command line. Unlike running from the binary, you will need to ensure you meet all the requirements outlined in the section below.
    2.  [ensemblator\_cli.py] - The CLI version of the *Ensemblator*. The installation is the same as for the GUI version, simply download, mark as executable, and ensure you meet all the requirements. The usage is different, in that this version operates from the command line, and thus is more amenable to being part of a automated pipeline.

### Requirements

#### Strict:

To use the *Ensemblator* source code, you will require Python 2 to be installed on your system. Furthermore, you will need the following python packages:

-   numpy
-   biopython
-   matplotlib
-   SciPy
-   scikit-learn

Some of these packages might be difficult to install using pip, but an alternative could be to use a scientific python installation like Anaconda.

#### Optional:

-   [muscle]

    > This software is needed for doing sequence alignments when building ensembles. This feature is VERY useful, I highly recommend it. Make sure that it is installed, and in your path as ‘muscle’. If you have muscle in your path, but are still encountering errors, please try running from the command line. Sometimes when clicking the icon from the desktop, the PATH variable does not get imported correctly. I don’t really know why this happens.

### Usage:

#### The main menu:

![image][1]

##### Prepare Input

Clicking this will open a window for preparing input for analysis.

  [image]: screenshots/all3.png
  [Installation]: #installation
  [Requirements]: #requirements
  [Usage]: #usage
  [The main menu]: #the-main-menu
  [Preparing your files for analysis]: #preparing-your-files-for-analysis
  [Analyzing your prepared ensemble]: #analyzing-your-prepared-ensemble
  [Understanding the output]: #understanding-the-output
  [Known Bugs]: #known-bugs
  [Bug Reports]: #bug-reports
  [ensemblator.py]: ensemblator.py
  [ensemblator\_cli.py]: ensemblator_cli.py
  [muscle]: http://www.drive5.com/muscle/
  [1]: screenshots/main_menu.png
