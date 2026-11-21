# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
try:
    from docutils import nodes
    nodes.reprunicode = lambda x: x # Workaround for bug in nbsphinx-link
except:
    pass # docutils may not be installed for test suit

new_path = ['..\\', '..\\thermosteam\\', '..\\Bioindustrial-Park\\', '..\\How2STEAM\\', '..\\hensmith\\']
for p in new_path:
     sys.path.insert(0, os.path.abspath(p))

# -- Project information -----------------------------------------------------

project = 'BioSTEAM'
copyright = '2018-2022, BioSTEAM Development Group'
author = 'Yoel Cortes-Pena'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''

# -- Images ------------------------------------------------------------------
# from sphinx.builders.html import StandaloneHTMLBuilder

# new_supported_image_types = [
#     'image/svg+xml',
#     'image/gif',
#     'image/png',
#     'image/jpeg'
# ]

# # construct it this way so that if Sphinx adds default support for additional images, such
# # as HEIC, then what we do is add any of those to the end. We start with the ones
# # we want to support in this order, then subtract them from the defaults to identify
# # any remaining items that we append to the end of the list

# additional_default_supported_images = list(set(StandaloneHTMLBuilder.supported_image_types) - set(new_supported_image_types))
# StandaloneHTMLBuilder.supported_image_types = new_supported_image_types + additional_default_supported_images

# -- General configuration ---------------------------------------------------

add_module_names = False

autodoc_member_order = 'bysource'

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx.ext.napoleon',
    'sphinxcontrib.jquery', # Some extensions require jquery, which used to come with Sphinx before v6.
    'matplotlib.sphinxext.plot_directive',
    'IPython.sphinxext.ipython_console_highlighting', # Fixes bug with pygments highlighting in nbsphinx (a workaround)
    'sphinx_design',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    # 'sphinx.ext.ifconfig',
    'sphinx.ext.autosummary',
    'sphinx_multitoc_numbering',
    'sphinx_autodoc_typehints',
    'myst_parser',
    'nbsphinx',
    'nbsphinx_link',
]

from functools import reduce
import biosteam as bst
import typing
import types
from types import UnionType, GenericAlias
from thermosteam import units_of_measure 
from numpy import typing as numpy_typing
from sphinx.util.typing import restify
defs = numpy_typing.__dict__ | units_of_measure.__dict__ | bst.__dict__ | typing.__dict__ | types.__dict__
defs['tmo'] = __import__('thermosteam')

def replace_forward_ref(annotation):
    if hasattr(annotation, '__forward_arg__'):
        forward_arg = annotation.__forward_arg__
        try: return eval(forward_arg, defs)
        except: return annotation
    elif isinstance(annotation, GenericAlias):
        return annotation
    elif isinstance(annotation, UnionType):
        return reduce(lambda a, b: a | b, [replace_forward_ref(i) for i in annotation.__args__])
    elif hasattr(annotation, '__args__'):
        try: annotation.__args__ = tuple([replace_forward_ref(i) for i in annotation.__args__])
        except: pass
        return annotation
    elif isinstance(annotation, str):
        try:
            return eval(annotation, defs)
        except:
            annotation = annotation.replace('(', '[').replace(')', ']')
            try:
                return eval(annotation, defs)
            except:
                return annotation
    else:
        return annotation

NoneType = type(None)

def typehints_formatter(annotation, config):
    # Called by sphinx_autodoc_typehints before its own rendering. Forward 
    # references are resolved using biosteam's namespace, optional types 
    # are rendered as "T, optional", and annotations are rendered with Sphinx 
    # (so nested annotations are not passed to this function again). 
    # Returning None falls back to the default rendering of 
    # sphinx_autodoc_typehints.
    resolved = replace_forward_ref(annotation)
    if isinstance(resolved, str): return None # Could not resolve
    mode = 'fully-qualified-except-typing' if config.typehints_fully_qualified else 'smart'
    suffix = ''
    if typing.get_origin(resolved) in (typing.Union, UnionType):
        args = typing.get_args(resolved)
        if NoneType in args:
            args = tuple([i for i in args if i is not NoneType])
            resolved = args[0] if len(args) == 1 else typing.Union[args]
            suffix = ', optional'
    try:
        return restify(resolved, mode) + suffix
    except Exception:
        return None

typehints_fully_qualified = False
always_use_bars_union  = True
simplify_optional_unions = False
always_document_param_types = False
typehints_document_rtype = False
typehints_use_rtype = False
autodoc_typehints_format = 'short'
imgmath_latex_preamble = r'\usepackage{xcolor}'
nbsphinx_execute = 'never'

# Specify the baseurls for the projects I want to link to
intersphinx_mapping = {
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'python': ('https://docs.python.org/3', None),
    'qsdsan': ('https://qsdsan.readthedocs.io/en/latest/', None),
    'thermo': ('https://thermo.readthedocs.io/', None),
    'chaospy': ('https://chaospy.readthedocs.io/en/master/', None),
}

# Allow exceptions to occur in notebooks
nbsphinx_allow_errors = True

# Do not show all members of the class
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ['.txt', '.rst', '.md']

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = [
    '_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints',
    # Not documentation pages
    'README.rst', 'requirements.txt', 'API/thermosteam/paper/**',
]

if os.environ.get("SPHINX_FAST_BUILD"):
    html_search_index = False
    # Skips all slow autodoc generation
    exclude_patterns.extend([
        "API/**",
        "tutorial/**",
        "**.ipynb",
        "**.py",
    ])

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.

html_theme = "pydata_sphinx_theme"

#
html_theme_options = {
    "logo" : {
        'image_light': 'logo.png',
        'image_dark': 'logo_dark.png'
    },
    "show_toc_level": 2,
    "announcement": (
        "<p>Workshop on Uncertainty, Sensitivity, and Optimization in March 2027, "
        "<a href='https://docs.google.com/spreadsheets/d/1QlP7iSZrLHkUxhk-2R8_QWomh3qHfkchJK1EvM5cpyo/edit?usp=sharing'>sign up here</a>!</p>"
    ),
    "external_links": [
      {"name": "Bioindustrial-Park", "url": "https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park"},
      {"name": "How2STEAM", "url": "https://mybinder.org/v2/gh/BioSTEAMDevelopmentGroup/How2STEAM/HEAD"},
      {"name": "QSDsan", "url": "https://qsdsan.readthedocs.io/en/latest/"},
      {"name": "YouTube Channel", "url": "https://www.youtube.com/@yoelcortes-pena2100/videos"},
  ]
}


html_sidebars = {
    "tutorial/index": [],
    "contributing/index": [],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["images", '_static']

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/custom.css',
]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'BioSTEAMdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'BioSTEAM.tex', 'BioSTEAM Documentation',
     'Yoel Cortes-Pena and Guest Group', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'biosteam', 'BioSTEAM Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'BioSTEAM', 'BioSTEAM Documentation',
     author, 'BioSTEAM', 'The Bioprocess Simulation and Techno-Economic Analysis Modules',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for markdown extentions -----------------------------------------

myst_heading_anchors = 3

# -- Design inputs -----------------------------------------------------------

# Design inputs of units are documented twice in the class docstring:
#
# * In the "Parameters" section, since they can be passed as key-word 
#   arguments (with an ", optional" suffix in the type). Parameters already 
#   documented have priority; otherwise, the description is given by 
#   `<design input>._make_doc(param=True)`.
#
# * In a "Design inputs" section at the end of the docstring (right before 
#   other attributes), as attributes documented by `<design input>.__doc__`:
#
#     Design inputs
#
#     T: float
#         Operating temperature [K]. T must be positive.
#
# Design inputs are not listed again as members of the class. 
import re
import inspect
from sphinx.util.typing import stringify_annotation
from biosteam._unit import DesignInput

def get_design_inputs(obj):
    if not isinstance(obj, type): return []
    names = getattr(obj, '_design_inputs', None)
    if not names: return []
    design_inputs = []
    for name in names:
        design_input = getattr(obj, name, None)
        if isinstance(design_input, DesignInput): design_inputs.append(design_input)
    return design_inputs

def format_type(design_input):
    return_type = getattr(design_input, '_return_type', None)
    if return_type is None: return None
    return stringify_annotation(return_type, 'smart')

def is_numpy_header(lines, i):
    return (
        i + 1 < len(lines) and lines[i].strip()
        and not lines[i].startswith((' ', '\t'))
        and re.fullmatch(r'-{3,}', lines[i + 1].strip()) is not None
    )

def add_design_input_parameters(app, what, name, obj, options, lines):
    # Runs before napoleon converts numpy sections to reST
    if what != 'class': return
    design_inputs = get_design_inputs(obj)
    if not design_inputs: return
    N = len(lines)
    header = next((i for i in range(N) if lines[i].strip() == 'Parameters' and is_numpy_header(lines, i)), None)
    if header is None:
        # Add a parameters section before the first section (or at the end)
        index = next((i for i in range(N) if is_numpy_header(lines, i)), N)
        lines[index:index] = ['Parameters', '----------', '']
        start = end = index + 2
    else:
        start = header + 2
        end = next((i for i in range(start, len(lines)) if is_numpy_header(lines, i)), len(lines))
    keys = {i.__name__: i for i in design_inputs}
    documented = set()
    for index in range(start, end):
        line = lines[index]
        if not line.strip() or line.startswith((' ', '\t')): continue
        head, colon, type = line.partition(':')
        names = [i.split('=', 1)[0].strip().lstrip('*') for i in head.split(',')]
        documented.update(names)
        if len(names) == 1 and names[0] in keys:
            # Design inputs are optional key-word arguments
            type = type.strip() or format_type(keys[names[0]]) or ''
            if not type.endswith('optional'): 
                type = type + ', optional' if type else 'optional'
            lines[index] = f'{head.rstrip()} : {type}'
    entries = []
    for design_input in design_inputs:
        key = design_input.__name__
        if key in documented: continue
        type = format_type(design_input)
        entries.append(f'{key} : {type}, optional' if type else f'{key} : optional')
        entries.append('    ' + design_input._make_doc(param=True))
    if not entries: return
    while end > start and not lines[end - 1].strip(): end -= 1 
    lines[end:end] = entries + ([''] if end < len(lines) else [])

def add_design_inputs_section(app, what, name, obj, options, lines):
    # Runs after napoleon (and sphinx_autodoc_typehints) converts the docstring to reST
    if what != 'class': return
    design_inputs = get_design_inputs(obj)
    if not design_inputs: return
    section = ['']
    for design_input in design_inputs:
        doc = design_input.__doc__ or ''
        section.append(f'.. py:attribute:: {design_input.__name__}')
        type = format_type(design_input)
        if type: section.append(f'   :type: {type}')
        section.append('')
        section.extend(['   ' + i.strip() if i.strip() else '' for i in doc.splitlines()])
        section.append('')
    # Add at the end of the docstring (after references, if any), 
    # right before the members of the class (e.g., attributes)
    while lines and not lines[-1].strip(): lines.pop()
    lines.extend(section)

# -- Thermodynamic models ----------------------------------------------------

# Properties listed in a "Thermodynamic models" section of a class docstring 
# (e.g., thermosteam.Chemical) are not documented again as members.
models_section = 'Thermodynamic models'
models_cache = {}

def get_owner(function):
    # Return the class that defines the function (if any)
    names = getattr(function, '__qualname__', '').split('.')[:-1]
    if not names or '<locals>' in names: return None
    obj = sys.modules.get(getattr(function, '__module__', None))
    for i in names: obj = getattr(obj, i, None)
    return obj if isinstance(obj, type) else None

def get_model_names(cls):
    if cls in models_cache: return models_cache[cls]
    lines = (cls.__dict__.get('__doc__') or '').splitlines()
    names = set()
    for i, line in enumerate(lines):
        if line.strip() == models_section and i + 1 < len(lines) and set(lines[i + 1].strip()) == {'-'}:
            indent = len(line) - len(line.lstrip())
            for line in lines[i + 2:]:
                if not line.strip(): continue
                current = len(line) - len(line.lstrip())
                if current < indent: break
                if current > indent: continue # Description
                if set(line.strip()) == {'-'}: break # Next section
                name = re.match(r'\s*(\w+)', line)
                if name: names.add(name.group(1))
            break
    models_cache[cls] = names
    return names

def skip_model_members(name, obj):
    function = getattr(obj, 'fget', obj)
    owner = get_owner(function)
    return owner is not None and name in get_model_names(owner)

def skip_design_input_members(app, what, name, obj, skip, options):
    # Design inputs are already documented in the "Design inputs" section
    if isinstance(obj, DesignInput): return True
    # Thermodynamic models are already documented in the class docstring
    if skip_model_members(name, obj): return True

def short_unit_signature(app, what, name, obj, options, signature, return_annotation):
    # Design inputs and other key-word only parameters of units are 
    # summarized as **kwargs (they are documented in the docstring)
    if what != 'class' or not (isinstance(obj, type) and issubclass(obj, bst.Unit)): return
    try:
        parameters = [*inspect.signature(obj).parameters.values()]
    except (TypeError, ValueError):
        return
    short = []
    kwargs = False
    for i in parameters:
        if i.kind in (i.KEYWORD_ONLY, i.VAR_KEYWORD): 
            kwargs = True
            continue
        short.append(i.replace(annotation=i.empty))
    if kwargs: short.append(inspect.Parameter('kwargs', inspect.Parameter.VAR_KEYWORD))
    return str(inspect.Signature(short)), return_annotation

def setup(app):
    app.connect('autodoc-process-signature', short_unit_signature, priority=100) # Before sphinx_autodoc_typehints (first result is used)
    app.connect('autodoc-process-docstring', add_design_input_parameters, priority=100)
    app.connect('autodoc-process-docstring', add_design_inputs_section, priority=900)
    app.connect('autodoc-skip-member', skip_design_input_members)
