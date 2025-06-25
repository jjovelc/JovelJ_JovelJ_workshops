# Getting Ready for the Workshop: Shotgun Metagenomics Data Analysis

This document provides basic instructions to help you prepare for the upcoming workshop. If you have any questions or run into issues, please feel free to contact the workshop instructors.

## Table of Contents

1. [Get an ARC Account](#get-an-arc-account)
2. [Install mamba in ARC](#install-mamba-in-arc)
3. [Install and Configure VSCode](install_config_VSCode.md)
4. [Setting Up R in VSCode](#setting-up-R-in-vscode)
5. [Basic Linux Tutorial (optional)](https://www.juanjovel.ca/blogs/HPC_tutorial.html)
6. [Basic R Tutorial (optional)](https://www.juanjovel.ca/blogs/R_tutorial.html)


---
## Get an ARC Account

Processing metagenomic libraries often requires access to a High-Performance Computing (HPC) environment such as **ARC**, the research computing cluster at the University of Calgary.

If you do not already have an ARC account, please ask your PI (Principal Investigator) to request one on your behalf. Below is a template that your PI can send to the ARC support team at: `support@hpc.ucalgary.ca`.

```
Subject: Request for ARC User Accounts

Dear ARC Support Team,

My name is XXXXX, and I am a faculty member in the Faculty of Veterinary Medicine at the University of Calgary.

As part of my research activities, I require access to the high-performance computing resources provided by ARC. I would like to request the creation of user accounts for the following individuals:

xxxxxxxxxx@ucalgary.ca
xxxxxxxxxx@ucalgary.ca

Thank you for your assistance and support. Please let me know if you need any further information.

Best regards,
```
---

## Install mamba in ARC

**Note:** If you already have Conda or Mamba installed, you can skip this step.

Mamba is a package management system similar to Conda, but faster. The primary difference is that Conda is written in Python, whereas Mamba is written in C++, making it more efficient for certain tasks.

### Steps for Installation:

1. **Download the installer** from the official Mamba repository:

[Download Mambaforge](https://github.com/conda-forge/miniforge/releases)

Alternatively, you can use `wget` to download the installer directly:


```bash
wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh

```

2. Run the installer:

```bash
bash Miniforge3-25.3.0-3-Linux-x86_64.sh
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

---

## Setting Up R in VSCode

After installing VSCode, following the instructions provided [here](install_config_VSCode.md), please proceed to install R and configure it in VSCode.

### 📦 1. Install R

🔹 On macOS (via Homebrew)

```bash
brew install --cask r
```

Or download from: https://cran.r-project.org/

🔹 On Windows

Download and install from:
👉 https://cran.r-project.org/bin/windows/base/

### 🧩 3. Install Required VS Code Extensions
Open VS Code, then:

Press Ctrl+Shift+X (or Cmd+Shift+X on macOS)

Search and install the following extensions:

  - R Language – by Yuki Ueda
  - (Optional) R Debugger – by Manuel Hentschel


### ▶️  4. Using R in VS Code

Start R Terminal
Press Ctrl+Shift+P → choose R: Create R Terminal

This will launch an R session in the VS Code terminal. Install the following package:

Run Code
Open an .R script

Select lines and press Ctrl+Enter (Windows) or Cmd+Enter (Mac) to send to the terminal

- Install the following packages (copy these commands in the R terminal):

install.packages("devtools")

install.packages('ape')

install.packages(c("tidyverse", "data.table", "httpgd", "languageserver"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

BiocManager::install("metagenomeSeq")

BiocManager::install("DESeq2")

BiocManager::install("lefser")

install.packages("remotes")
remotes::install_github("vegandevs/vegan")





