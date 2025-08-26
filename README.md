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

## Source code to edge detection experiment

**Prerequisites**: The KITT (Kermit Image Toolkit) collection is needed for this project to work. The folders needed are loaded in **"setup.m"** file, where paths can be configured. Please download the files from [KITT repository](https://github.com/giaracvi/KITT).

To execute the experiment, run the file **"superLauncher.m"** where the OS (win, mac, linux) must be chosen in order to build the correct path syntax.

The parameters configuration is located in the file **"infoMaker.m"**. The source and data paths must also be configured according to your folders location.

Each one of the phases of the experiment is located in one file, as follows:

- **"aioMaker.m"** contains all of the instructions of the exection of the experiment.
  
- **"cpCollecter.m"** takes individual statistical results and collects them all in order to have a global result of the dataset

- **"README.md"** this file.

--------------------------------------------------------------------------------
## Source code to classification problem
This code is based on FARC-HD algorithm, public available at KEEL dataset repository (https://github.com/SCI2SUGR/KEEL).

To execute the expirment, build a Java project. The configuration files are available in the Sugeno folder. Precisely, the Param.txt will define the folds to be used in training and test. We provide the well-known iris dataset in folder data. In the param.txt file, consider:
Type of Inference = (WR (0), AC(1) or Sugeno/Choquet (2), Probabilistic Sum (3));
Type of fuzzy measure = The Fuzzy mesure used with the operator (Cardinality (1), Dirac (2), Wmean (3), OWA (4) and PM (5));
Type of aggregation = The FG-functional used in the FRM.

Finally, we highlight tha the FRM is available in the file RuleBase.java. You may find it in FARC-HD folder.
