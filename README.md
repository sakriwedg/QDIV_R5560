# QDIV Data Acquisition

This repository contains the software for controlling the data acquisition of **Charge Division (QDIV) detectors** using the **CAEN R5560** module.

## Installation

Before running the acquisition software, create a Python virtual environment and install the required dependencies.

### Windows

Open a terminal in the repository directory and run:

```bash
python -m venv .venv
.venv\Scripts\activate
```

### macOS / Linux

Open a terminal in the repository directory and run:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### Install the requirements

Once the virtual environment is activated, install the required Python packages with:

```bash
pip install -r requirements.txt
```

## Running the software

With the virtual environment activated, start the QDIV acquisition GUI with:

```bash
python QDIV_gui.py
```

The GUI provides the interface for controlling the data acquisition through the CAEN R5560 module.

