# Quantification of scRNAseq libraries

## Installation of mamba.

Mamba is a package management system that is analogous to Conda, just faster. The reason is, Conda is writen in Python, while Mamba is writen in C++. 

1. Download installer from:

https://github.com/conda-forge/miniforge/releases

```bash
wget https://github.com/conda-forge/miniforge/releases/download/24.7.1-2/Mambaforge-24.7.1-2-Linux-x86_64.sh
```

2. Run the executable to install Mamba:

```bash
bash Mambaforge-24.7.1-2-Linux-x86_64.sh
```

Follow instructions and allow mamba to modify your bashrc (yes in last question).

3. Type command 'mamba' if the command is not available, type que following command:

```bash
echo 'export PATH=/work/vetmed_data/mamba/bin/:$PATH' >> ~/.bashrc 
```

4. source ~/.bashrc:

```bash
source ~/.bashrc
```


## Installation of Kallisto-Bus


