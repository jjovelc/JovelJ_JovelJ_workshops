# VSCode for HPC: Installation and Remote Connection Tutorial

## Table of Contents
1. [Installing VSCode](#installing-vscode)
2. [Essential Extensions for HPC Work](#essential-extensions-for-hpc-work)
3. [Setting Up SSH Keys](#setting-up-ssh-keys)
4. [Connecting to HPC via SSH](#connecting-to-hpc-via-ssh)
5. [File Transfer Methods](#file-transfer-methods)
6. [Working with Remote Files](#working-with-remote-files)
7. [Terminal Usage](#terminal-usage)
8. [Troubleshooting Common Issues](#troubleshooting-common-issues)

## Installing VSCode

### Windows Installation

1. **Download VSCode**
   - Go to [https://code.visualstudio.com/](https://code.visualstudio.com/)
   - Click "Download for Windows"
   - Choose the appropriate version (User Installer recommended)

2. **Run the Installer**
   - Double-click the downloaded `.exe` file
   - Accept the license agreement
   - Choose installation location (default is fine)
   - **Important**: Check these boxes during installation:
     - "Add to PATH" (enables command line usage)
     - "Create a desktop icon"
     - "Add 'Open with Code' action to Windows Explorer file context menu"

3. **Complete Installation**
   - Click "Install" and wait for completion
   - Launch VSCode when installation finishes

### Mac Installation

1. **Download VSCode**
   - Go to [https://code.visualstudio.com/](https://code.visualstudio.com/)
   - Click "Download for Mac"
   - Choose Universal or your specific chip (Apple Silicon/Intel)

2. **Install VSCode**
   - Open the downloaded `.zip` file
   - Drag "Visual Studio Code" to your Applications folder
   - Launch VSCode from Applications or Spotlight

3. **Add to PATH (Optional but Recommended)**
   - Open VSCode
   - Press `Cmd+Shift+P` to open Command Palette
   - Type "shell command" and select "Shell Command: Install 'code' command in PATH"
   - This allows you to open VSCode from terminal with `code .`

## Essential Extensions for HPC Work

Install these extensions for optimal HPC experience:

### Required Extensions

1. **Remote - SSH**
   - Extension ID: `ms-vscode-remote.remote-ssh`
   - Enables SSH connections to remote servers
   - Install: Go to Extensions (Ctrl+Shift+X), search "Remote SSH"

2. **Remote - SSH: Editing Configuration Files**
   - Extension ID: `ms-vscode-remote.remote-ssh-edit`
   - Helps edit SSH configuration files

### Recommended Extensions

3. **Remote Explorer**
   - Usually comes with Remote SSH
   - Provides GUI for managing remote connections

4. **Python** (if you work with Python)
   - Extension ID: `ms-python.python`
   - Full Python support with IntelliSense

5. **Jupyter** (for notebook support)
   - Extension ID: `ms-toolsai.jupyter`
   - Run Jupyter notebooks remotely

6. **GitLens** (for Git integration)
   - Extension ID: `eamodio.gitlens`
   - Enhanced Git capabilities

## Setting Up SSH Keys

SSH keys provide secure, password-free authentication to your HPC system.

### Generate SSH Key Pair

**Windows (using Command Prompt or PowerShell):**
```bash
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
```

**Mac (using Terminal):**
```bash
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
```

**During key generation:**
- Press Enter for default file location (`~/.ssh/id_rsa`)
- Enter a secure passphrase (optional but recommended)
- Press Enter to confirm passphrase

### Copy Public Key to HPC

**Method 1: Using ssh-copy-id (Mac/Linux)**
```bash
ssh-copy-id username@hpc-server.edu
```

**Method 2: Manual Copy (Windows/Mac)**
```bash
# Display your public key
cat ~/.ssh/id_rsa.pub

# Copy the output, then SSH to your HPC and add it to authorized_keys
ssh username@hpc-server.edu
mkdir -p ~/.ssh
echo "your_public_key_here" >> ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
chmod 700 ~/.ssh
```

## Connecting to HPC via SSH

### Configure SSH Connection

1. **Open VSCode**
2. **Access Remote Explorer**
   - Click the Remote Explorer icon in the sidebar (computer with plug icon)
   - Or press `Ctrl+Shift+P` and type "Remote SSH: Connect to Host"

3. **Add New SSH Target**
   - Click the "+" icon next to "SSH TARGETS"
   - Enter connection string: `ssh username@hpc-server.edu`
   - Choose SSH config file location (usually `~/.ssh/config`)

### Edit SSH Configuration

For more advanced setups, edit your SSH config file:

1. **Open SSH Config**
   - Press `Ctrl+Shift+P`
   - Type "Remote SSH: Open SSH Configuration File"
   - Select your config file

2. **Add HPC Configuration**
```
Host myhpc
    HostName hpc-cluster.university.edu
    User your_username
    Port 22
    IdentityFile ~/.ssh/id_rsa
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

3. **Connect to HPC**
   - In Remote Explorer, click the folder icon next to your HPC host
   - Or press `Ctrl+Shift+P` and select "Remote SSH: Connect to Host"
   - Choose your HPC from the list
   - VSCode will open a new window connected to the HPC

## File Transfer Methods

### Method 1: Using VSCode File Explorer

**Once connected to HPC:**
1. **Navigate Remote Files**
   - Use the Explorer panel to browse HPC file system
   - Right-click folders/files for context menu options

2. **Upload Files**
   - Right-click in Explorer panel
   - Select "Upload" to upload files from local computer
   - Or drag and drop files from local file explorer

3. **Download Files**
   - Right-click on remote files
   - Select "Download" to save to local computer

### Method 2: Using Integrated Terminal

**Open terminal in VSCode:**
- Press `Ctrl+`` (backtick) or Terminal → New Terminal

**Transfer commands:**
```bash
# Upload file to HPC (run from local terminal)
scp local_file.txt username@hpc-server.edu:~/remote_directory/

# Upload directory (recursive)
scp -r local_directory/ username@hpc-server.edu:~/remote_directory/

# Download file from HPC (run from local terminal)
scp username@hpc-server.edu:~/remote_file.txt ./local_directory/

# Download directory
scp -r username@hpc-server.edu:~/remote_directory/ ./local_directory/
```

**Using rsync for large transfers:**
```bash
# Sync local to remote
rsync -avz --progress local_directory/ username@hpc-server.edu:~/remote_directory/

# Sync remote to local
rsync -avz --progress username@hpc-server.edu:~/remote_directory/ ./local_directory/
```

### Method 3: Git Integration

If your code is in Git repositories:
```bash
# Clone repository on HPC
git clone https://github.com/username/repository.git

# Push changes from local
git add .
git commit -m "Update code"
git push

# Pull changes on HPC
git pull origin main
```

## Working with Remote Files

### Opening Remote Directories

1. **Connect to HPC** (as described above)
2. **Open Folder**
   - File → Open Folder (Ctrl+K Ctrl+O)
   - Navigate to your desired directory on HPC
   - Click "OK"

### Editing Remote Files

- Files open directly in VSCode editor
- All changes are saved directly to HPC
- Use Ctrl+S to save as usual
- Full IntelliSense and extension support available

### Managing Multiple Connections

- Each remote connection opens in a new VSCode window
- You can have multiple HPC connections simultaneously
- Local VSCode instance remains separate from remote instances

## Terminal Usage

### Accessing HPC Terminal

**In remote VSCode window:**
1. **Open Terminal**
   - Terminal → New Terminal (Ctrl+`)
   - Terminal opens directly on HPC system

2. **Multiple Terminals**
   - Click "+" to create additional terminal sessions
   - Use dropdown to switch between terminals

### Common HPC Commands

```bash
# Check job queue
squeue -u $USER

# Submit job
sbatch job_script.sh

# Check system resources
sinfo

# Load modules
module avail
module load python/3.9

# Check disk usage
du -sh *
df -h

# Monitor running processes
htop
```

### Running Interactive Sessions

```bash
# Request interactive node
srun --pty --time=01:00:00 --mem=8G bash

# Or using interactive queue
salloc --time=01:00:00 --mem=8G
```

## Troubleshooting Common Issues

### Connection Problems

**Issue: "Could not establish connection"**
- Verify HPC server address and username
- Check if SSH keys are properly configured
- Test connection in regular terminal: `ssh username@hpc-server.edu`
- Check if HPC allows SSH connections from your network

**Issue: "Permission denied (publickey)"**
- Ensure public key is in `~/.ssh/authorized_keys` on HPC
- Check key permissions: `chmod 600 ~/.ssh/id_rsa`
- Verify SSH agent is running: `ssh-add -l`

### Performance Issues

**Slow file operations:**
- Use compression for large transfers: `scp -C`
- Consider rsync for incremental transfers
- Work with smaller files when possible

**VSCode server crashes:**
- Check HPC disk quota: `quota -u $USER`
- Clear VSCode server cache: `rm -rf ~/.vscode-server`
- Restart connection

### Extension Issues

**Extensions not working remotely:**
- Some extensions need to be installed on remote server
- Go to Extensions panel and look for "Install in SSH: hostname"
- Python/Jupyter extensions may need remote installation

### File Permissions

**Permission denied errors:**
```bash
# Fix common permission issues
chmod 700 ~/.ssh
chmod 600 ~/.ssh/authorized_keys
chmod 600 ~/.ssh/id_rsa
chmod 644 ~/.ssh/id_rsa.pub
```

## Best Practices

### Security
- Always use SSH keys instead of passwords
- Use strong passphrases for SSH keys
- Regularly update your SSH keys
- Don't share private keys

### Performance
- Close unused remote connections
- Use workspace settings for project-specific configurations
- Enable auto-save for remote files (File → Auto Save)
- Use .gitignore to exclude large data files from version control

### Workflow
- Keep local and remote code synchronized via Git
- Use relative paths in your code for portability
- Create shell scripts for common HPC tasks
- Document your setup for team members

## Advanced Tips

### SSH Config Optimization
```
Host myhpc
    HostName hpc-cluster.university.edu
    User username
    Port 22
    IdentityFile ~/.ssh/id_rsa
    ControlMaster auto
    ControlPath ~/.ssh/control_%h_%p_%r
    ControlPersist 1h
    ServerAliveInterval 60
    ServerAliveCountMax 3
    Compression yes
```

### Port Forwarding for Jupyter
```bash
# Forward Jupyter port from HPC to local machine
ssh -L 8888:localhost:8888 username@hpc-server.edu
```

### Using VSCode Tasks
Create `.vscode/tasks.json` for common operations:
```json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "Submit Job",
            "type": "shell",
            "command": "sbatch",
            "args": ["${file}"],
            "group": "build"
        }
    ]
}
```

This tutorial should get you started with using VSCode as your primary interface for HPC development and data management. The combination of remote editing, integrated terminal, and seamless file transfer makes VSCode an excellent choice for HPC workflows.
