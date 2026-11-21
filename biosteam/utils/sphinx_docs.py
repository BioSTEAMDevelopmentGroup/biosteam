# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Render text (or the documentation of python objects) with Sphinx using the
configuration of BioSTEAM's documentation (docs/conf.py).
"""
import os
import sys
import inspect
import tempfile
from io import StringIO

__all__ = ('sphinx_docs',)

#: Default configuration directory (docs folder of the biosteam repository).
default_confdir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'docs')
)

def autodoc_directive(obj, members=True, inherited_members=False):
    if members:
        options = '   :members:\n'
        if inherited_members: options += '   :inherited-members:\n'
    else:
        options = ''
    if inspect.isclass(obj):
        if issubclass(obj, BaseException): directive = 'autoexception'
        else: directive = 'autoclass'
    elif inspect.ismodule(obj):
        return f'.. automodule:: {obj.__name__}\n' + options
    elif inspect.isroutine(obj):
        directive = 'autofunction'
    else:
        raise TypeError(
            f'cannot generate documentation for {type(obj).__name__!r} objects'
        )
    module = obj.__module__
    qualname = obj.__qualname__
    if '.' in qualname and directive == 'autofunction': directive = 'automethod'
    text = f'.. {directive}:: {module}.{qualname}\n'
    if directive in ('autoclass', 'autoexception'): text += options
    return text

def sphinx_docs(source, builder='text', format='rst', confdir=None,
                intersphinx=False, warnings=None, confoverrides=None,
                members=True, inherited_members=False):
    """
    Return documentation rendered by Sphinx as text using the settings in
    BioSTEAM's docs/conf.py (extensions, napoleon, type hints, design inputs,
    etc.).

    Parameters
    ----------
    source : str or object
        Text to render (reStructuredText or Markdown) or a python object
        (class, function, method, or module) to document with autodoc.
    builder : str, optional
        Sphinx builder (e.g., 'text', 'html'). Defaults to 'text'.
    format : str, optional
        Format of the text, either 'rst' or 'md'. Defaults to 'rst'.
        Ignored if `source` is a python object.
    confdir : str, optional
        Directory with the conf.py file. Defaults to the docs folder of the
        biosteam repository.
    intersphinx : bool, optional
        Whether to load intersphinx inventories (requires internet).
        Defaults to False.
    warnings : file-like, optional
        Stream for Sphinx warnings. Defaults to sys.stderr.
    confoverrides : dict, optional
        Configuration values to override.
    members : bool, optional
        Whether to document the members (e.g., attributes, properties,
        design inputs, and methods) of classes and modules. Defaults to True.
        Ignored if `source` is a string.
    inherited_members : bool, optional
        Whether to also document inherited members. Defaults to False.
        Ignored if `source` is a string or `members` is False.

    Examples
    --------
    >>> import biosteam as bst
    >>> print(bst.utils.sphinx_docs('**Hello** :math:`x^2`')) # doctest: +SKIP
    **Hello** x^2
    >>> print(bst.utils.sphinx_docs(bst.Flash)) # doctest: +SKIP

    """
    from sphinx.application import Sphinx
    from sphinx.util.docutils import docutils_namespace, patch_docutils
    if confdir is None: confdir = default_confdir
    confdir = os.path.abspath(confdir)
    if not os.path.isfile(os.path.join(confdir, 'conf.py')):
        raise FileNotFoundError(
            f'no conf.py file in {confdir!r}; the biosteam docs folder is only '
            'available when biosteam is installed from the repository '
            '(pass `confdir` to use another directory)'
        )
    if isinstance(source, str):
        text = source
        if format not in ('rst', 'md'):
            raise ValueError(f"format must be either 'rst' or 'md', not {format!r}")
    else:
        text = autodoc_directive(source, members, inherited_members)
        format = 'rst'
    overrides = {
        'exclude_patterns': [],
        'suppress_warnings': ['toc.not_included', 'config.cache'],
    }
    if not intersphinx: overrides['intersphinx_mapping'] = {}
    if confoverrides: overrides.update(confoverrides)
    if warnings is None: warnings = sys.stderr
    with tempfile.TemporaryDirectory() as tmp:
        srcdir = os.path.join(tmp, 'src')
        outdir = os.path.join(tmp, 'out')
        doctreedir = os.path.join(tmp, 'doctrees')
        os.mkdir(srcdir)
        with open(os.path.join(srcdir, 'index.' + format), 'w', encoding='utf-8') as file:
            file.write(text)
        cwd = os.getcwd()
        try:
            # Isolate docutils registries (as in sphinx-build) so that
            # Sphinx can be run multiple times in the same process
            with patch_docutils(confdir), docutils_namespace():
                app = Sphinx(
                    srcdir, confdir, outdir, doctreedir, builder,
                    confoverrides=overrides, status=StringIO(), warning=warnings,
                    freshenv=True,
                )
                app.build(force_all=True)
        finally:
            os.chdir(cwd)
        filename = os.path.join(outdir, 'index' + app.builder.out_suffix)
        with open(filename, encoding='utf-8') as file:
            return file.read()
