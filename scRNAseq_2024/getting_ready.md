# Getting ready for the scRNAseq data analysis course

This document guides you through the installation of all software required for the scRNAseq data analysis course

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
echo 'export PATH=/work/vetmed_data/mamba/bin/:$PATH' >> ~/.bashrc 
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



## Installation of Kallisto-Bus (in ARC server)

### 1. Create a new mamba environment to host kallisto and bus tools:

```bash
mamba create -n kallisto_bus python=3.12 -y
``` 

### 2. Install kallisto

```bash
# Activate environment
mamba activate kallisto_bus

# Install kallisto
mamba install -c bioconda kallisto -y
```

Try kallisto. If properly installed, you should see the following after typing 'kallisto':

---
```
kallisto 0.51.1

Usage: kallisto `<CMD>` [arguments] ...

Where <CMD> can be one of:

    index         Builds a kallisto index 
    quant         Runs the quantification algorithm 
    quant-tcc     Runs quantification on transcript-compatibility counts
    bus           Generate BUS files for single-cell data 
    h5dump        Converts HDF5-formatted results to plaintext
    inspect       Inspects and gives information about an index
    version       Prints version information
    cite          Prints citation information

    Running kallisto <CMD> without arguments prints usage information for <CMD>
```
---

### 3. Install kb-python tools (in ARC server)

```bash
pip install kb-python
```

Try kb-python. If properly installed you should see the following after typing kb:

---
```
usage: kb [-h] [--list] <CMD> ...

kb_python 0.28.2

positional arguments:
  <CMD>
    info      Display package and citation information
    compile   Compile kallisto and bustools binaries from source
    ref       Build a kallisto index and transcript-to-gene mapping
    count     Generate count matrices from a set of single-cell FASTQ files

options:
  -h, --help  Show this help message and exit
  --list      Display list of supported single-cell technologies
```
---

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

### 4. Install DoubletFinder

```bash

install.packages("DoubletFinder")

```

### 5. Install EnhacedVolcano

```bash

install.packages("EnhancedVolcano")

```
