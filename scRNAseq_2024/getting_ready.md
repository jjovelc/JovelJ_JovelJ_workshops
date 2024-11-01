# Getting ready for the scRNAseq data analysis course

This document guides you through the installation of all software required for the scRNAseq data analysis course. If you are using a Mac laptop, you are ready to go. If you are using a Windows laptop, you would need to install a terminal, options are MobaXterm [install MobaXterm] (https://mobaxterm.mobatek.net/download.html), Windows subsystem for linux [install WSL] (https://learn.microsoft.com/en-us/windows/wsl/install), or PuTTY [install PuTTY] (https://www.putty.org/). If you are using a Linux laptop, you are a crack and don't need any assistance.

## Installation of Mamba (in ARC server)

**Note:** If you already have Conda or Mamba installed, you can skip this step.

Mamba is a package management system similar to Conda, but faster. The primary difference is that Conda is written in Python, whereas Mamba is written in C++, making it more efficient for certain tasks.

### Steps for Installation:

1. **Download the installer** from the official Mamba repository:

   [Download Mambaforge](https://github.com/conda-forge/miniforge/releases)

   Alternatively, you can use `wget` to download the installer directly:


```bash
wget https://github.com/conda-forge/miniforge/releases/download/24.7.1-2/Mambaforge-24.7.1-2-Linux-x86_64.sh

```

2. Run the installer:

```bash
bash Mambaforge-24.7.1-2-Linux-x86_64.sh
```

Follow prompts and allow mamba to modify your .bashrc file (yes in last question).

3. Type command 'mamba' if the command is not available, type que following command:

```bash
echo 'export PATH=/home/<your-user-name>/mambarforge:$PATH' >> ~/.bashrc 
```

4. source .bashrc:

```bash
source ~/.bashrc
```

5. Initialize mamba:

```bash
mamba init
```

6. Source .bashrc again

```bash
source ~/.bashrc
```

After completing these steps, Mamba should be installed and ready to use. You can verify the installation by typing mamba in the terminal.


## Installation of salmon

```bash
mamba create -n salmon -c bioconda salmon -y
``` 

If correctly installed, you can activate the environment and test salmon:

```bash
mamba activate salmon

salmon

salmon v1.10.3

'Usage:  salmon -h|--help or 
        salmon -v|--version or 
        salmon -c|--cite or 
        salmon [--no-version-check] <COMMAND> [-h | options]

'Commands:
     'index      : create a salmon index
     'quant      : quantify a sample
     'alevin     : single cell analysis
     'swim       : perform super-secret operation
     'quantmerge : merge multiple quantifications into a single file
```

You can also inspect the help of the alevin module, which is the one we will use for quantificaiton of our scRNAseq libraries.


```bash
salmon alevin -h

Version Server Response: Not Found

alevin
==========
salmon-based processing of single-cell RNA-seq data.

'alevin options:


'mapping input options:
  -l [ --libType ] arg                  Format string describing the library 
                                        type
  -i [ --index ] arg                    salmon index
  -r [ --unmatedReads ] arg             List of files containing unmated reads 
                                        of (e.g. single-end reads)
  -1 [ --mates1 ] arg                   File containing the #1 mates
  -2 [ --mates2 ] arg                   File containing the #2 mates


'alevin-specific Options:
  -v [ --version ]                      print version string
  -h [ --help ]                         produce help message
  -o [ --output ] arg                   Output quantification directory...
```

## Installation of R packages (in your laptop)

### 0. Install R and Rstudio

R and Rstudio installers can be downloaded from [R & Rstudio] (https://posit.co/download/rstudio-desktop/)


### 1. Install remotes

Run these commands from the R console. After installation has completed without errors, you can test a package by simply running library(package-name). IF the library was properly installed, no error message will be deployed, otherwise, the system would report a message like "Error in library(package-name) : there is no package called ‘package-name’". 

``` bash

install.packages("remotes")

# Load library (if installation successful)
library(remotes)

```

### 2. Install Seurat v5

```bash

remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

```

### 3. Install Azimuth

```bash

remotes::install_github('satijalab/azimuth', ref = 'master')

```

### 3. Install tidyverse

```bash

install.packages("tidyverse")

```

NOTE: If you run into any trouble, just contact me (juan.jovel@ucalgary.ca).


### 5. Install EnhacedVolcano

```bash

install.packages("EnhancedVolcano")

```

### 6. Install SeuratDisk

```bash

remotes::install_github("mojaveazure/seurat-disk")

```
