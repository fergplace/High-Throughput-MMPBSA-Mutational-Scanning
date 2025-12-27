# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
from pathlib import Path
import os
sys.path.append(Path(os.path.join(Path(os.getcwd()).parent, "HTMS_Amber")).as_posix())


project_root = Path(__file__).resolve().parent.parent.parent 
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
    print(f"Sphinx DEBUG: Added project root '{project_root}' to sys.path.")
else:
    print(f"Sphinx DEBUG: Project root '{project_root}' already in sys.path.")

autodoc_mock_imports = ["modeller"]

project = 'HTMS in Amber'
copyright = '2025, Fergus Place'
author = 'Fergus Place'
release = '0.1.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # For Google/NumPy style docstrings
    'sphinx.ext.viewcode', # To link to source code
    'sphinx.ext.autosummary', # For generating summary tables
    'myst_parser',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
]
#mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
#new for md as landing
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
master_doc = 'index'
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
]
mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"], ["$", "$"]],
        "displayMath": [["\\[", "\\]"], ["$$", "$$"]],
        "processEscapes": True,
        "packages": {'[+]': ['ams', 'color', 'physics', 'cancel', 'extpfeil', 'mhchem', 'ams', 'braket', 'unicode']}, # You can add more packages here if needed for specific LaTeX commands
    },
}

templates_path = ['_templates']
#exclude_patterns = ['old_index.rst']
exclude_patterns = ['old_index.md', "old_index copy.rst"]

# html_theme_options = {
#     # 1. GitHub Icon Link in the top-right navbar
#     "icon_links": [
#         {
#             "name": "GitHub",
#             "url": "https://github.com/fergplace/HTMS_Amber",
#             "icon": "fa-brands fa-square-github",
#             "type": "fontawesome",
#         }
#     ],
    
#     # 2. Navbar & Logo
#     "logo": {
#         "text": "HTMS in Amber",
#     },
#     "navbar_align": "left",
    
#     # 3. Dark Mode / Light Mode toggle
#     # PyData theme includes a toggle button by default. 
#     # The styles below set the code highlighting for each mode.
#     "pygments_light_style": "tango",
#     "pygments_dark_style": "monokai",
    
#     # 4. Footer/Social
#     "footer_start": ["copyright"],
#     "footer_end": ["sphinx-version", "theme-version"],
# }



html_theme_options = {
    # 1. GitHub Icon Link in the top-right navbar
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/fergplace/HTMS_Amber",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        }
    ],
    
    # 2. Navbar & Header Layout
    "logo": {
        "text": "HTMS in Amber",
    },
    "navbar_align": "left",
    # CRITICAL: Setting this to 0 moves all page links out of the top header 
    # and into the left sidebar for a cleaner look.
    "header_links_before_dropdown": 2, 
    
    # 3. Sidebar and Content Navigation
    # Control the left (primary) and right (secondary) sidebars
    "show_nav_level": 1,
    "navigation_depth": 4,
    "collapse_navigation": False,
    "show_toc_level": 2, # Show sections of the current page in the right sidebar
    
    # 4. Dark Mode / Light Mode toggle
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
    
    "secondary_sidebar_items": ["page-toc"],
    # 5. Footer/Social
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version", "theme-version"],
}

# Define what appears in the sidebars
# html_sidebars = {
#     "**": [
#         "sidebar-nav-bs",     # Global site navigation (Left)
#     ]
# }

html_sidebars = {
    # For the API section (assuming your modules/API use .rst)
    #"modules/**": ["sidebar-nav-bs", "sourcelink"],
    "HTMS_Amber/**": ["sidebar-nav-bs", "sourcelink"],
    
    # For everything else (User Guides/MD), only show navigation
    "**": ["sidebar-nav-bs"]
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_baseurl = 'https://fergplace.github.io/HTMS_Amber/'
#html_theme = 'alabaster' #old one
html_theme = "pydata_sphinx_theme" #testing it 
html_static_path = ['_static']
# html_css_files = [
#     'custom.css',
# ]
autosummary_generate = True