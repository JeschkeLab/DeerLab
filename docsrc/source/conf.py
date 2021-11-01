# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------

import sys
import os

sys.path.append(os.path.abspath('.'))

# Import Read-the-Docs (RTD) theme
from sphinx.locale import _
from sphinx_rtd_theme import __version__

# Project details
project = 'DeerLab'
copyright = '2019-2021, Luis Fábregas-Ibáñez, Stefan Stoll, and others'
author = 'Fabregas Ibanez'
language = 'en'

# Print the HTML code for the landing page with dynamically compiled version number
version = open(os.path.join('..','..','VERSION')).read().splitlines()[0]
string = f"""
:raw-html:`<div class="topfloatcontainer">
        <div class="illustration">
            <img src="_static/landingpage.svg" alt="DeerLab Illustration">
        </div>
        <div class="title">
            <p class="titlep"><span class="title1">DeerLab</span><span class="title2">Docs</span><span class="version">{version}</span></p>
            <span style="color:#586069; font-size:20px; line-height:1.5;">
            DeerLab is a comprehensive free scientific software package for Python focused on modelling, penalized least-squares regression, and uncertainty quantification. It also provides specialized models and tools for the analysis of dipolar EPR (electron paramagnetic resonance) spectroscopy techniques such as DEER (double electron-electron resonance), and others. 
            </span>
        <br>
        <a href="./user_guide.html" class="btn-quick", style="margin-top:20px;">   Get started → </a>
        <br>
        </div>
    </div>`
"""
import re
string = re.sub('\s+',' ',string)
rst_epilog = f""" 
.. role:: raw-html(raw)
   :format: html

.. |title_version| replace:: {string}
"""

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# Add sphinx extensions
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
    'm2r2',
    'sphinx_copybutton'
]

#sys.path.append('../../deerlab')
#autosummary_mock_imports = ['deerlab']
add_module_names = False
autosummary_generate = True  # Turn on sphinx.ext.autosummary

from sphinx_gallery.sorting import ExplicitOrder
sphinx_gallery_conf = {
    'filename_pattern': 'ex_',
    'examples_dirs': '../../examples',   # path to your example scripts
    'gallery_dirs': 'auto_examples',  # path to where to save gallery generated output
    'remove_config_comments': True,
    'download_all_examples': False,
    'subsection_order': ExplicitOrder(['../../examples/basic',
                                       '../../examples/advanced',
                                       '../../examples/other']),
}

copybutton_prompt_text = ">>> "

# Warnings suppression
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
exclude_patterns = ['.', './functions']
numpydoc_show_class_members = False
# Render Latex math equations as svg instead of rendering with JavaScript
imgmath_image_format = 'png' if os.name=='nt' else 'svg'
imgmath_dvisvgm = 'dvisvgm'
imgmath_latex_preamble = r'''
\newcommand{\mr}[1]{\mathrm{#1}}
\newcommand{\mx}[1]{\boldsymbol{#1}}
\newcommand{\vc}[1]{\boldsymbol{#1}}
\DeclareMathOperator*{\argmin}{\arg\!\min}
'''

# Setup template stuff
templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'
intersphinx_mapping = {
    'sphinx': ('http://www.sphinx-doc.org/en/stable/', None),
}
#html_theme = 'sphinx_rtd_theme'
html_theme = "pydata_sphinx_theme"

# Integrate version control system
# -------------------------------------------------------------
html_context = {
    "display_github": False, # Integrate GitHub
    "github_user": "JeschkeLab", # Username
    "github_repo": "DeerLab", # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/source/", # Path in the checkout to the docs root
    'version' : version,                                  
}


# Read-the-Docs options configuration
# -------------------------------------------------------------
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/JeschkeLab/DeerLab",
            "icon": "fab fa-github",
        },
        {
            "name": "Twitter",
            "url": "https://twitter.com",
            "icon": "fab fa-twitter",
        },
        {
            "name": "Issues",
            "url": "https://github.com/JeschkeLab/DeerLab/issues",
            "icon": "fas fa-bug",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/DeerLab/",
            "icon": "fas fa-cube",
        },
    ],
    #"navbar_end": ["navbar-icon-links.html", "search-field.html"],
    #"navbar_align": "left"
}
#html_sidebars = {
#  "**": []
#}
html_sidebars = {
    "index": []
}
html_copy_source = False
html_theme_path = ["../.."]
html_show_sourcelink = True
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']

# Extensions to theme docs
def setup(app):
    from sphinx.domains.python import PyField
    from sphinx.util.docfields import Field
    #app.add_css_file('/source/_static/custom.css')
    #app.add_stylesheet('/source/_static/theme_override.css')
    app.add_object_type(
        'confval',
        'confval',
        objname='configuration value',
        indextemplate='pair: %s; configuration value',
        doc_field_types=[
            PyField(
                'type',
                label=_('Type'),
                has_arg=False,
                names=('type',),
                bodyrolename='class'
            ),
            Field(
                'default',
                label=_('Default'),
                has_arg=False,
                names=('default',),
            ),
        ]
    )
    # Patch the MATLAB lexer to correct wrong highlighting
    from sphinx.highlighting import lexers
    from pygments.lexers import python
    from pygments.token import Name
    from pygments.filters import NameHighlightFilter
    from deerlab import dd_models, bg_models
    import matplotlib.pyplot as plt 
    import deerlab
    import numpy
    import inspect
    #python_lexer = python.PythonLexer()
    #dd_functions = [o[0] for o in inspect.getmembers(bg_models) if inspect.isfunction(o[1])]
    #bg_functions = [o[0] for o in inspect.getmembers(dd_models) if inspect.isfunction(o[1])]
    #pl_functions = [o[0] for o in inspect.getmembers(plt) if inspect.isfunction(o[1])]
    #dl_functions = [o[0] for o in inspect.getmembers(deerlab) if inspect.isfunction(o[1])]
    #np_functions = [o[0] for o in inspect.getmembers(numpy) if inspect.isfunction(o[1])]

    #python_lexer.add_filter(NameHighlightFilter(
    #        names= dl_functions + dd_functions + bg_functions + pl_functions + np_functions,
    #        tokentype=Name.Function,
    #        ))
    #app.add_lexer('python', python_lexer)

# These folders are copied to the documentation's HTML output
html_static_path = ['_static']

# Add path to custom CSS file to overwrite some of the default CSS settings
html_css_files = [
#    'custom.css',
#    'table_util.css',
#    'table_main.css',
    'theme_override.css'
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Default role
default_role = 'math'  # with this, :math:`\psi` can be written simply as `\psi`


# -- Options for HTML output -------------------------------------------------
#html_theme = "sphinx_rtd_theme"
#html_theme_path = ["_themes", ]
html_static_path = ['_static']
html_title = 'DeerLab'
highlight_language = 'python'
primary_domain = 'py'
html_logo = '_static/logo_docs.svg'


# Design pygments patch for MATLAB true code highlighting
# --------------------------------------------------------

# Import the pygments python library
from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, Number, Operator, Generic, Text, Other, Comment, Whitespace

class MyFancyStyle(Style):
    default_style = "default"
    styles = {
        Text:                   '#fff',
        Comment:                '#6a737d',
        Keyword:                '#d73a49',
        Operator.Word:          '#d73a49',
        Name.Variable:          '#dddee4',
        Name.Function:          '#9065e2',
        Name.Class:             '#0000FF',
        Name.Builtin:           '#9065e2',
        Name.Attribute:         '#9065e2',
        String:                 '#244679',
        Operator:               '#005cc5',
        Number:                 '#0e56b1',
    }

# Create patch for applying the style 	
def pygments_monkeypatch_style(mod_name, cls):
    import sys
    import pygments.styles
    cls_name = cls.__name__
    mod = type(__import__("os"))(mod_name)
    setattr(mod, cls_name, cls)
    setattr(pygments.styles, mod_name, mod)
    sys.modules["pygments.styles." + mod_name] = mod
    from pygments.styles import STYLE_MAP
    STYLE_MAP[mod_name] = mod_name + "::" + cls_name

pygments_monkeypatch_style("my_fancy_style", MyFancyStyle)
pygments_style = "my_fancy_style"