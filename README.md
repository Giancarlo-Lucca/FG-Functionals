# FG-Functionals
This repository is related with the paper (a, b)-FG-functionals: a generalization of the Sugeno integral with floating domains in arbitrary closed real intervals and its applications submited to the journal 

This paper is autored by:
- Tiago Asmus,
- Giancarlo Lucca,
- Cedric Marco-Detchart,
- Joaquín Arellano,
- Humberto Bustince,
- Graçaliz Pereira Dimuro

We provide two different folders each one related with a different application:

1 - Experimental approach to edge detection

2 - Experimental approach to classification problem

--------------------------------------------------------------------------------

## Source code for edge detection experiment

**Prerequisites**: The KITT (Kermit Image Toolkit) collection is needed for this project to work. The folders needed are loaded in **"setup.m"** file, where paths can be configured. Please download the files from [KITT repository](https://github.com/giaracvi/KITT).

To execute the experiment, run the file **"superLauncher.m"** where the OS (win, mac, linux) must be chosen in order to build the correct path syntax.

The parameters configuration is located in the file **"infoMaker.m"**. The source and data paths must also be configured according to your folders location.

Each one of the phases of the experiment is located in one file, as follows:

- **"aioMaker.m"** contains all of the instructions of the exection of the experiment.
  
- **"cpCollecter.m"** takes individual statistical results and collects them all in order to have a global result of the dataset

- **"README.md"** this file.
