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

new_path = ['..\\', '..\\thermosteam\\', '..\\Bioindustrial-Park\\', '..\\How2STEAM\\']
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
    'sphinxcontrib.jquery', # Some extensions require jquery, which used to come with Sphinx before v6.
    "matplotlib.sphinxext.plot_directive",
    'IPython.sphinxext.ipython_console_highlighting', # Fixes bug with pygments highlighting in nbsphinx (a workaround)
    'sphinx_design',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
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

try:
    import sphinx_autodoc_typehints as sat
    def format_annotation(annotation, config) -> str:  # noqa: C901 # too complex
        # Special cases
        if isinstance(annotation, sat.ForwardRef):
            return annotation.__forward_arg__
        if annotation is None or annotation is type(None):  # noqa: E721
            return ":py:obj:`None`"
        if annotation is Ellipsis:
            return ":py:data:`...<Ellipsis>`"
    
        if isinstance(annotation, tuple):
            return sat.format_internal_tuple(annotation, config)
    
        try:
            module = sat.get_annotation_module(annotation)
            class_name = sat.get_annotation_class_name(annotation, module)
            args = sat.get_annotation_args(annotation, module, class_name)
        except ValueError:
            return str(annotation).strip("'")
    
        # Redirect all typing_extensions types to the stdlib typing module
        if module == "typing_extensions":
            module = "typing"
    
        full_name = f"{module}.{class_name}" if module != "builtins" else class_name
        fully_qualified: bool = getattr(config, "typehints_fully_qualified", False)
        prefix = "" if fully_qualified or full_name == class_name else "~"
        if module == "typing" and class_name in sat._PYDATA_ANNOTATIONS:
            role = "data"
        else:
            role = "class"
        args_format = "\\[{}]"
        formatted_args: str | None = ""
    
        # Some types require special handling
        if full_name == "typing.NewType":
            args_format = f"\\(``{annotation.__name__}``, {{}})"
            role = "class" if sys.version_info >= (3, 10) else "func"
        elif full_name == "typing.TypeVar":
            params = {k: getattr(annotation, f"__{k}__") for k in ("bound", "covariant", "contravariant")}
            params = {k: v for k, v in params.items() if v}
            if "bound" in params:
                params["bound"] = f" {format_annotation(params['bound'], config)}"
            args_format = f"\\(``{annotation.__name__}``{', {}' if args else ''}"
            if params:
                args_format += "".join(f", {k}={v}" for k, v in params.items())
            args_format += ")"
            formatted_args = None if args else args_format
        elif full_name == "typing.Optional":
            args = tuple(x for x in args if x is not type(None))  # noqa: E721
        elif full_name in ("typing.Union", "types.UnionType") and type(None) in args:
            if len(args) == 2:
                full_name = "typing.Optional"
                args = tuple(x for x in args if x is not type(None))  # noqa: E721
            else:
                simplify_optional_unions: bool = getattr(config, "simplify_optional_unions", True)
                if not simplify_optional_unions:
                    full_name = "typing.Optional"
                    args_format = f"\\[:py:data:`{prefix}typing.Union`\\[{{}}]]"
                    args = tuple(x for x in args if x is not type(None))  # noqa: E721
        elif full_name == "typing.Callable" and args and args[0] is not ...:
            fmt = [format_annotation(arg, config) for arg in args]
            formatted_args = f"\\[\\[{', '.join(fmt[:-1])}], {fmt[-1]}]"
        elif full_name == "typing.Literal":
            formatted_args = f"\\[{', '.join(repr(arg) for arg in args)}]"
        elif full_name == "types.UnionType":
            return " | ".join([format_annotation(arg, config) for arg in args])
    
        if args and not formatted_args:
            try:
                iter(args)
            except TypeError:
                fmt = [format_annotation(args, config)]
            else:
                fmt = [format_annotation(arg, config) for arg in args]
            formatted_args = args_format.format(", ".join(fmt))
    
        result = f":py:{role}:`{prefix}{full_name}`{formatted_args}"
        return result
    
    def typehints_formatter(annotation, config):
        # original = str(annotation).replace("'", "")
        name = format_annotation(annotation, config)
        if name.startswith('Optional['):
            name = name.replace('Optional[', '')[:-1] + ', optional'
        elif name.startswith(':py:data:`~typing.Optional`\['):
            name = name.replace(':py:data:`~typing.Optional`\[', '')[:-1] + ', optional'
        if not name.startswith(':py:class:') and not name.startswith(':py:data:'):
            pyclass_added = False
            for i in ['Iterable', 'Sequence', 'Collection', 'str', 'list', 'tuple', 'dict', 'int',
                      'bool', 'set', 'frozenset', 'float']:
                name = name.replace(i, f':py:class:`{i}`')
                pyclass_added = True
            for i in ['Unit', 'Facility', 'HeatUtility', 'PowerUtility', 'System', 'TEA', 'HXutility']:
                name = name.replace(i, f':py:class:`~biosteam.{i}`')
                pyclass_added = True
            for i in ['Chemical', 'Thermo', 'Stream', 'ReactionSet']:
                name = name.replace(i, f':py:class:`~thermosteam.{i}`')
                pyclass_added = True
            if pyclass_added: 
                name = name.replace('[', '\[').replace(':py:class::py:class:', ':py:class:')
        name = name.replace('bst.', '').replace('tmo.', '').replace(
            '`biosteam', '`~biosteam').replace('._unit', '')
        # file = os.path.join(os.path.dirname(__file__), 'annotations.txt')
        # with open(file, 'a') as f: f.write(f"{name}  ({original})\n")
        return name
except:
    pass

simplify_optional_unions = True
always_document_param_types = False
typehints_document_rtype = False
typehints_use_rtype = False
autodoc_typehints_format = 'short'
imgmath_latex_preamble = r'\usepackage{xcolor}'
nbsphinx_execute = 'never'

# Specify the baseurls for the projects I want to link to
intersphinx_mapping = {
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'python': ('https://docs.python.org/3', None),
    'qsdsan': ('https://qsdsan.readthedocs.io/en/latest/', None),
    'thermo': ('https://thermo.readthedocs.io/', None),
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
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

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
    # "announcement": (
    #     "<p>Join us on Feb 17, 9:15-10:15am CST, for a BioSTEAM workshop! "
    #     "<a href='mailto: biosteamdevelopmentgroup@gmail.com'>Email us for details</a></p>"
    # ),
    "external_links": [
      {"name": "Bioindustrial-Park", "url": "https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park"},
      {"name": "How2STEAM", "url": "https://mybinder.org/v2/gh/BioSTEAMDevelopmentGroup/How2STEAM/HEAD"},
      {"name": "QSDsan", "url": "https://qsdsan.readthedocs.io/en/latest/"},
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





