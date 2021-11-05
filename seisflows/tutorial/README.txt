SEISFLOWS3 DEVEL VERSION INVERSION EXAMPLE README FOR MAUI
STARTED: 27/08/2021
LAST UPDATED: 5/11/2021
CONTACT: Bryant Chow (bhchow@alaska.edu)

################################################################################

INTRODUCTION

################################################################################

This README is an introduction to the SeisFlows3 (seisflows) package. 
It discusses the command line interface (CLI), the structure of a seisflows
working directory, and runs a basic inversion.

The steps are in logical order and most provide some commands to run to set up
the example. Some are optional, which means they only provide information but 
contain no critical commands to run.

################################################################################

STEP 1: Environment setup

################################################################################

(1) SeisFlows needs to be run in a Conda environment. This must be done on the
Maui ancillary node. To get there you need to SSH from maui using:

    $ ssh -X w-mauivlab01.maui.nesi.org.nz

(2) We will use my pre-built Conda environment to run SeisFlows3. 

    $ source activate /nesi/project/gns03247/PyPackages/conda_envs/seisflows

NOTE: You can also create your own conda environment using the text file 
provided which lists all the packages in my environment. This is not necessary
    
    $ cd /scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/seisflows/examples/TUTORIAL/misc/conda_envs
    $ conda create --name seisflows seisflows_conda_env.txt

(3) You should now be able to access the CLI by typing 'seisflows' in the terminal
which will list out the available commands in the Seisflows package

    $ seisflows

usage: seisflows [-h] [-w [WORKDIR]] [-p [PARAMETER_FILE]]
                 [--path_file [PATH_FILE]]
                 {setup,modules,configure,init,submit,resume,restart,clean,par,check,convert,reset,inspect,debug,edit}
                 ...

================================================================================

                     SeisFlows: Waveform Inversion Package

================================================================================

optional arguments:
  -h, --help            show this help message and exit
  -w [WORKDIR], --workdir [WORKDIR]
                        The SeisFlows working directory, default: cwd
  -p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Parameters file, default: 'parameters.yaml'
  --path_file [PATH_FILE]
                        Legacy path file, default: 'paths.py'

command:
  Available SeisFlows arguments and their intended usages

    setup               Setup working directory from scratch
    modules             Print available module names
    configure           Fill parameter file with defaults
    init                Initiate working environment
    submit              Submit initial workflow to system
    resume              Re-submit previous workflow to system
    restart             Remove current environment and submit new workflow
    clean               Remove active working environment
    par                 View and edit parameter file
    check               Check state of an active environment
    convert             Convert model file format
    reset               Reset underlying machinery
    inspect             View inheritenace and ownership
    debug               Start interactive debug environment
    edit                Open source code file in text editor

'seisflows [command] -h' for more detailed descriptions of each command.


################################################################################

STEP 2: Directory setup

################################################################################

(1) We can use the CLI to create a new SF3 working directory. You can choose 
your own (empty) directory to follow along. 

   $ cd /scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/seisflows/examples/
   $ seisflows setup

(2) A file named 'parameters.yaml' has been created. This is an ASCII par file
that controls everything within SF3 and will be your main point of control over
how SF3 runs inversions, submits jobs, and generally interacts with the cluster

(3) If you open 'parameters.yaml' you will see a few modules which define each
of the main components of SF3. To view the available module choices defined
by SF3 we can run the modules command

    $ seisflows modules

(4) Because we're working on Maui and interested in performing an inversion 
(or just only the first part of an inversion), we need to change some of the 
default module choices in the 'parameters.yaml' file. Using a text editor, 
set the module choices to:

WORKFLOW: inversion
SOLVER: specfem3d
SYSTEM: maui
OPTIMIZE: LBFGS
PREPROCESS: pyatoa
POSTPROCESS: base

(4) Now we can run the 'configure' command which will automatically fill in the
remainder of the parameter file based on the choices we've set for our modules

    $ seisflows configure

(5) Open up the 'parameters.yaml' file again and you will see that the configure
command has set up the entire SeisFlows3 parameter file with docstrings
explaining what each of the parameters does. 

Note: Parameters that need to be filled in by the user are denoted with
 '!!! REQUIRED PARAMETER !!!'
For a general user, these parameters will need to be filled in based on their 
specific cluster, problem, and intended outcome. 

(6) To save some time, I have pre filled a parameters.yaml file, you can see 
which values have been set to what by using vimdiff

    $ vimdiff parameters.yaml TUTORIAL/parameters_complete.yaml

(7) Let's use this pre-filled parameter file to save time

    $  cp -f TUTORIAL/parameters.yaml .

Great, the main, critical component of SeisFlows3, the parameters.yaml file is 
now complete. There are two more components that we require to run SeisFlows3,
the DATA/ directory, and the MODEL/ directory.
    

################################################################################

STEP 3: DATA AND SPECFEM3D WORKING DIRECTORY (OPTIONAL)

################################################################################

* This step explains how to setup the DATA directory for SeisFlows3,
which tells SF3 and SPECFEM3D what sources and stations need to be included
in the inversion. 

* When running earthquake forward simulations, SPECFEM3D requires a DATA/ 
directory which contains (among other things):

    (a) CMTSOLUTION
    (b) STATIONS
    (c) Par_file

* Similarly, SeisFlows3 requires a DATA/ directory with the same files.
These files are distributed by SeisFlows3 to each of the SPECFEM3D working 
directories under the hood.

* One key difference is that SeisFlows3 requires N number of CMTSOLUTION files,
where N is the number of earthquakes you want to simulate in parallel. N should
also match the parameters.yaml value: NTASK

You can check NTASK with the following command

    $ seisflows par NTASK

* Each CMTSOLUTION file should be different, and have a suffix to distinguish
themselves from one another. One example would be to use earthquake IDs from
the organization you collect your earthquakes from, another example would be 
numerical ordering (e.g., CMTSOLUTION_001, CMTSOLUTION_002, ..., CMTSOLUTION_N)

! NOTE: The following steps are not necessary for the tutorial, they are just 
used to describe the SPECFEM3D directory structure

(1) For this example problem, the data and initial model have already been 
generated to keep things simple. We'll have a look to see the different parts.

    $ ls -l TUTORIAL/specfem3d

(2) Here we have a slimmed down version of a SPECFEM3D working directory, which
only requires a file for the compiled binaries (bin/) and the DATA/ directory
    
    $ ls -l TUTORIAL/specfem3d/DATA

! NOTE: Whatever is in this DATA directory (defined as PATHS.SPECFEM_DATA in 
the parameters.yaml file) will get copied into all the SeisFlows working 
directories (scales by number of events). So if you have large files, e.g., 
*.xyz files from an external tomography model that are not necessary for the
inversion, you may want to create a separate DATA/ directory that doesn't 
contain these files to avoid unncessary storage requirements.

(3) Here we see the DATA/ directory has parts (a)--(c) as mentioned above, as 
well as the directory meshfem3D_files

    $ ls -l TUTORIAL/specfem3d/DATA/meshfem3D_files

(4) These mesh files define the structured grid and velocity model used for
this example. The mesh defines the North Island of New Zealand split up into
N=40 chunks. The velocity model is a simple 1D layered halfspace, based
on the velocity model in Ristau (2008).

(5) We have already generated the DATABASES_MPI/ files which define the velocity
model in terms of the SPECMFEM3D binary (.bin) files. You can see those here:

    
    $ ls -l TUTORIAL/specfem3d/OUTPUT_FILES/DATABASES_MPI

(6) We have also generated a target model which is identical to the initial 
model except all the velocity (vp, vs) values are reduced by 30%. We will
be using SeisFlows to attempt to recover the velocity changes between
the initial and final models

* So now that we have our initial model and data, we can start our inversion


################################################################################

STEP 4: INITATE SEISFLOWS3 (OPTIONAL)

################################################################################

* We first want to initiate a SeisFlows3 working directory and have a look at 
all the different directories

* To initiate SeisFlows we can use the init command, to see what it does we
can run with the -h (help) flag (i.e., $ sf init -h)

(1) The init command essentially runs the SeisFlows workflow up to the point
where the master job is submitted, allowing the user to see the working 
directory without interacting with the cluster. Let's run it:

    $ sf init
    $ ls output

(2) You can see that the 'init' command created the output/ directory which 
contains .p (pickle) files, and .json (ascii) files. These files represent
the state of the SeisFlows workflow: all of tits initiated modules, and the
paths and parameters that are set by the User.

    $ cat output/seisflows_parameters.json

(3) If we open the parameters file we can see all the information that was 
listed in our parmameters.yaml file, which are accessed by SeisFlows to
control how it operates


################################################################################

STEP 5: DEBUG MODE (OPTIONAL)

################################################################################

* One of the most important tools in SeisFlows is its debug mode. This 
interactive iPython environment lets you explore an active workflow, run
individual functions/commands, and debug errors.

(1) Activate the debug mode, and when prompted type 'n' to get to iPython

    $ sf debug
    > n

(2) We are now in the usual IPython environment, except that the entire 
SeisFlows architecture is loaded into memory. We can explore the entire system,
but for now just as an example we'll just have a look at paths, parameters
and loaded modules

    > print(PAR)   # parameters are stored in a dict object
    > print(PATH)  # paths are also stored as dicts

* We can also access individual SeisFlows modules by calling their generic name

    In [6]: workflow
    Out[6]: <seisflows3.workflow.inversion.Inversion at 0x2aab3bb20550>

(3) We will use debug mode more in later tutorials. Any changes made here are
temporary unless saved using the 'workflow.checkpoint()' command, or if a 
function which writes to memory is called. To exit just type exit() at each
iPython prompt

    In [7]: exit()
    > exit()

(4) Just so we have a blank slate to run our inversion, we will run the 'clean'
command to get rid of all the created SeisFlows components

    $ sf clean


################################################################################

STEP 6: RUN THE INVERSION

################################################################################

* SeisFlows is an automated tool, which means that once we have set it up 
(parameters.yaml), we can just hit run and watch the output.

(1) Now we can run SeisFlows with the 'submit' command. When prompted, type 'y'
to submit the master job

    $ sf submit
    $ y

(2) Now that we have submitted the master job, we can check that it is running 
on the cluster

    $ squeue --account=gns03247

(3) We can also track the progress of the inversion by watching the output
logs. NOTE: ctrl+c to quit the tail command

    $ tail -f *log

* The output-%a.log file is the stdout where SeisFlows will tell you where it
is in the workflow. These correspond to various modules and functions within
the source code

(4) The inversion will likely take some time, in the meantime we can explore
the directory that the workflow has created


################################################################################

STEP 7: UNDERSTANDING THE SEISFLOWS WORKING DIRECTORY

################################################################################

* Once your inversion starts running, it will create a directory structure
to keep track of the workflow and all the data generated. It should look
something like this:

(seisflows) [ancil] chowbr@w-mauivlab01 [examples] $ ls
error-14907303.log  output/              output.logs/  output.stats/     scratch/
logs.old/           output-14907303.log  output.optim  parameters.yaml@  TUTORIAL/


(a) output-*.log
    ASCII text file where SeisFlows will print out status updates as it step 
    through the workflow. The number is the cluster job number

(b) error-*.log
    ASCII text file where SeisFlows will print traceback information if
    a workflow exits unexpectedly. If your job crashes, look here to figure out
    why

(c) output.optim
    ASCII text file where SeisFlows will print out information about the
    optimization algorithm, including iteration, step count and misfit info.

(d) output.logs/
    When individual processes are run (e.g., a simulation), they produce their
    own output files which are written here. Files are named after job numbers

(e) logs.old/
    If you re-run a workflow (e.g., using $ sf resume), new log files (a and b)
    will be created. The old ones are moved here so that the information
    is retained

(f) output.stats/
    List-like stats information related to the optimization algorithm are
    written out here, to be used for e.g., plotting or analysis

(g) output/
    Location where all permament results are stored. These include the pickle
    files and parameter files which store the state of the workflow (see 6.3)
    Similarly gradients, updated velocity models, waveforms etc. are stored 
    here if requested by the User.

(h) scratch/
    The location of all the large, temporary working files required to run the 
    workflow. Each of the SeisFlows modules has its own directory in scratch/

    (1) 
