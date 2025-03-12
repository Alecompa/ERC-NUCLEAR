# ERC-NUCLEAR Analysis Code
Shared analysis code and data for NUCLEAR project

## Requirements
- python3
- numpy
- scipy
- matplotlib
- nbstripout

## Installation
Clone the repository and install the analysis library in your python environment with:
```bash
pip install -e .  
```
## Usage
Simply import the analysis library in your notebook or script with
```python
import modules.analysis as analysis
```
Examples of usage script and notebooks for the different functions are given in ```src/notebooks``` and ```src/scripts``` folders.

## Setting up nbstripout for clean notebook commits

To prevent committing unnecessary output in Jupyter notebooks, we use [`nbstripout`](https://github.com/kynan/nbstripout). Please install and activate it before contributing:

### **Installation**
1. Install `nbstripout`:
   ```bash
   pip install nbstripout
   ```
2. Enable it for the repository:
   ```bash
   nbstripout --install --attributes .gitattributes
   ```
This will automatically strip output from .ipynb files before committing, keeping our repository clean.