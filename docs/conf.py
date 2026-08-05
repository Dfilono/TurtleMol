'''Sphinx configuration for TurtleMol.'''

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

project = 'TurtleMol'
copyright = '2023, Dominick Filonowich'
author = 'Dominick Filonowich'

try:
    from TurtleMol import __version__
except ImportError:
    __version__ = 'development'

release = __version__
version = __version__
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']
autodoc_member_order = 'bysource'
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_theme = 'alabaster'
html_title = 'TurtleMol'
html_logo = 'images/logo.png'
html_static_path = ['_static']
