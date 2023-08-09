# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------

import sys
import os
import re
from sphinx.locale import _
import warnings

# Get DeerLab version
version = open(os.path.join('..','..','VERSION')).read().splitlines()[0]

# Warnings suppression
warnings.filterwarnings("ignore", category=FutureWarning)

# Add path
sys.path.append(os.path.abspath('.'))

# Project details
project = 'DeerLab'
copyright = '2019-2023, Luis Fábregas-Ibáñez, Stefan Stoll, Hugo Karas and others'
author = 'Fabregas Ibanez'
language = 'en'

# Sphinx extensions
# ----------------------------------------------------------------------
extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinxcontrib.httpdomain',
    'sphinxcontrib.ghcontributors',
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.imgmath',
    'numpydoc',
    'sphinx_gallery.gen_gallery',
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx_issues',
    'sphinx_copybutton',
    'sphinx_design'
]

# Configureation of Sphinx-Issues
# ----------------------------------------------------------------------
# GitHub repo
issues_github_path = "JeschkeLab/deerlab"
# equivalent to
issues_uri = "https://github.com/JeschkeLab/deerlab/issues/{issue}"
issues_pr_uri = "https://github.com/JeschkeLab/deerlab/pull/{pr}"
issues_commit_uri = "https://github.com/JeschkeLab/deerlab/commit/{commit}"


# Configuration of Sphinx-Autosymmary
# ----------------------------------------------------------------------
add_module_names = False
# Turn on sphinx.ext.autosummary
autosummary_generate = True 
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# Configuration of the Sphinx-Gallery
# ----------------------------------------------------------------------
from sphinx_gallery.sorting import ExplicitOrder
sphinx_gallery_conf = {
    'filename_pattern': 'ex_(?!long_)',
    'examples_dirs': '../../examples',   # path to your example scripts
    'gallery_dirs': 'auto_examples',  # path to where to save gallery generated output
    'remove_config_comments': True,
    'download_all_examples': False,
    'subsection_order': ExplicitOrder(['../../examples/basic_simulations',
                                       '../../examples/basic',
                                       '../../examples/intermediate',
                                       '../../examples/advanced',
                                       '../../examples/other']),
    'capture_repr': ()
}

# Configuration of Sphinx-CopyButton
# ----------------------------------------------------------------------
copybutton_prompt_text = ">>> "



# Configuration of Latex-to-SVG
# ----------------------------------------------------------------------
# Render Latex math equations as svg instead of rendering with JavaScript
imgmath_image_format = 'png' if os.name=='nt' else 'svg'
imgmath_dvisvgm = 'dvisvgm'
imgmath_latex_preamble = r'''
\newcommand{\mr}[1]{\mathrm{#1}}
\newcommand{\mx}[1]{\boldsymbol{#1}}
\newcommand{\vc}[1]{\boldsymbol{#1}}
\DeclareMathOperator*{\argmin}{\arg\!\min}
'''

# Configuration of the HTML Theme Template
# ----------------------------------------------------------------------
templates_path = ['_templates']

# Setup template stuff
exclude_patterns = ['.', './functions']
numpydoc_show_class_members = False
html_theme = "furo"
source_suffix = '.rst'
exclude_patterns = []
master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'
intersphinx_mapping = {'sphinx': ('http://www.sphinx-doc.org/en/stable/', None)}
html_context = {
    "display_github": False, # Integrate GitHub
    "github_user": "JeschkeLab", # Username
    "github_repo": "DeerLab", # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/source/", # Path in the checkout to the docs root
    'version' : version,     
    "default_mode": "light",                             
}


html_copy_source = False
html_theme_path = ["../.."]
html_show_sourcelink = True
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']
html_theme_options = {
    "sidebar_hide_name": True,
    "light_logo": "logo_docs.svg",
    "dark_logo": "logo_docs_light.svg",
}


default_role = 'math'  # with this, :math:`\psi` can be written simply as `\psi`
html_title = 'DeerLab'
highlight_language = 'python'
primary_domain = 'py'

