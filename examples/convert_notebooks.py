import subprocess
import os

# def convert_notebooks_to_markdown(notebook_files):
#     for notebook in notebook_files:
#         if not os.path.exists(notebook):
#             print(f"Error: {notebook} not found.")
#             continue
            
#         print(f"Converting {notebook} to Markdown...")
#         try:
#             result = subprocess.run(
#                 ["jupyter", "nbconvert", "--to", "markdown", notebook],
#                 check=True,
#                 capture_output=True,
#                 text=True
#             )
#             print(f"Successfully converted: {notebook}")
            
#         except subprocess.CalledProcessError as e:
#             print(f"Failed to convert {notebook}.")
#             print(f"Error output: {e.stderr}")
#         except FileNotFoundError:
#             print("Error: 'jupyter' command not found. Please ensure nbconvert is installed.")
#             break

def convert_notebooks_to_markdown(notebook_files):
    """
    Converts a list of .ipynb files to .md format using jupyter nbconvert
    and places the output in the ../doc_maker/source directory.
    """
    # 1. Define the target directory (one level up, then into doc_maker/source)
    # We use abspath to ensure it calculates the parent correctly from the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(os.path.dirname(script_dir), "docs_maker", "source")
    
    for notebook in notebook_files:
        if not os.path.exists(notebook):
            print(f"Error: {notebook} not found in current folder.")
            continue
            
        print(f"Converting {notebook}...")
        
        try:
            # The --output-dir flag is the key to sending the .md file 
            # and its image folders to the parent directory path.
            subprocess.run(
                [
                    "jupyter", 
                    "nbconvert", 
                    "--to", "markdown", 
                    notebook, 
                    "--output-dir", target_dir
                ],
                check=True,
                capture_output=True,
                text=True
            )
            print(f"Success! File saved to: {target_dir}")
            
        except subprocess.CalledProcessError as e:
            print(f"Failed to convert {notebook}.")
            print(f"Error: {e.stderr}")
        except FileNotFoundError:
            print("Error: 'jupyter' not found. Run 'pip install nbconvert'.")
            break

if __name__ == "__main__":
    files_to_convert = [
        "how_to_install.ipynb",
        "how_to_run.ipynb",
        "Method_Overview.ipynb"
    ]
    
    convert_notebooks_to_markdown(files_to_convert)