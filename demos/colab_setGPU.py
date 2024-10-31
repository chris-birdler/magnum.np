import json
import argparse

def set_gpu_runtime(notebook_path):
    with open(notebook_path, "r") as f:
        notebook = json.load(f)
    
    # Add Colab GPU runtime metadata
    notebook['metadata']['colab'] = notebook.get('metadata', {}).get('colab', {})
    notebook['metadata']['colab']['hardwareAccelerator'] = 'GPU'
    
    # Save the notebook with updated metadata
    with open(notebook_path, "w") as f:
        json.dump(notebook, f, indent=2)

# Set up argument parsing
parser = argparse.ArgumentParser(description="Set Colab notebooks to use GPU runtime.")
parser.add_argument("notebook_paths", type=str, nargs='+', help="Paths to Jupyter notebooks (.ipynb files)")

# Parse arguments
args = parser.parse_args()

# Apply GPU runtime setting to each notebook
for notebook_path in args.notebook_paths:
    set_gpu_runtime(notebook_path)
