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
copyright = '2019-2021, Luis Fábregas-Ibáñez, Stefan Stoll, and others'
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
    'sphinx_copybutton'
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
    'filename_pattern': 'ex_',
    'examples_dirs': '../../examples',   # path to your example scripts
    'gallery_dirs': 'auto_examples',  # path to where to save gallery generated output
    'remove_config_comments': True,
    'download_all_examples': False,
    'subsection_order': ExplicitOrder(['../../examples/basic',
                                       '../../examples/advanced',
                                       '../../examples/other']),
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
# Setup template stuff
exclude_patterns = ['.', './functions']
numpydoc_show_class_members = False
html_theme = "pydata_sphinx_theme"
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
}
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/JeschkeLab/DeerLab",
            "icon": "fab fa-github",
        },
        {
            "name": "Discussions",
            "url": "https://github.com/JeschkeLab/DeerLab/discussions",
            "icon": "fas fa-comments",
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
    ]
}
html_sidebars = {
    "index": [],
    "modelsref": ["search-field"],
    "reference": ["search-field"],
    "_autosummary/**": ["search-field"],
    "examples": ["search-field"],
    "auto_examples/**": ["search-field"],
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

# These folders are copied to the documentation's HTML output
html_static_path = ['_static']

# Add path to custom CSS file to overwrite some of the default CSS settings
html_css_files = [
    'theme_override.css'
]
# Default role
default_role = 'math'  # with this, :math:`\psi` can be written simply as `\psi`
html_title = 'DeerLab'
highlight_language = 'python'
primary_domain = 'py'
html_logo = '_static/logo_docs.svg'


# Patch Code highlighting
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


# Print the HTML code for the landing page with dynamically compiled version number
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
string = re.sub('\s+',' ',string)
rst_epilog = f""" 
.. role:: raw-html(raw)
   :format: html

.. |title_version| replace:: {string}

.. |fix| replace:: :raw-html:`<span class="badge changelog_fix">Fix</span>`

.. |efficiency| replace:: :raw-html:`<span class="badge changelog_efficiency">Efficiency</span>`

.. |enhancement| replace:: :raw-html:`<span class="badge changelog_enhancement">Enhancement</span>`

.. |feature| replace:: :raw-html:`<span class="badge changelog_feature">Feature</span>`

.. |api| replace:: :raw-html:`<span class="badge changelog_api">API Change</span>`
"""
