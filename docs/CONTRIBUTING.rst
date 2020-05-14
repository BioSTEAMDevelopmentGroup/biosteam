Contributing to BioSTEAM
========================

General Process
---------------

Here’s the short summary, complete TOC links are below:

#. If you are a first-time contributor:

   * Go to https://github.com/yoelcortes/biosteam and click the “fork” button to create your own copy of the project.

   * Clone the project to your local computer::
    
        git clone https://github.com/your-username/biosteam.git
    
   * Change the directory::
    
        cd biosteam
    
   * Add the upstream repository::
    
        git remote add upstream https://github.com/yoelcortes/biosteam.git
    
   * Now, git remote -v will show two remote repositories named "upstream" (which refers to the biosteam repository), and "origin" (which refers to your personal fork).

#. Develop your contribution:

   * Pull the latest changes from upstream::

       git checkout master
       git pull upstream master

   * Create a branch for the feature you want to work on. Since the branch name will appear in the merge message, use a sensible name such as "Pump-design-enhancement"::

       git checkout -b Pump-design-enhancement

   * Commit locally as you progress (git add and git commit) Use a properly formatted commit message, write tests that fail before your change and pass afterward, run all the tests locally. Be sure to document any changed behavior in docstrings, keeping to the NumPy docstring standard.

#. To submit your contribution:

   * Push your changes back to your fork on GitHub::

       git push origin Pump-design-enhancement

   * Enter your GitHub username and password (repeat contributors or advanced users can remove this step by connecting to GitHub with SSH).

   * Go to GitHub. The new branch will show up with a green Pull Request button. Make sure the title and message are clear, concise, and self- explanatory. Then click the button to submit it.

   * If your commit introduces a new feature or changes functionality, post in https://github.com/yoelcortes/biosteam/issues to explain your changes. For bug fixes, documentation updates, etc., this is generally not necessary, though if you do not get any reaction, do feel free to ask for review.

Testing
-------

For a full assessment of BioSTEAM, you can run the following:

.. code-block:: python

   >>> from biosteam.tests import test_biosteam
   >>> test_biosteam()

`test_biosteam` runs all the `doctests <https://docs.python.org/3.6/library/doctest.html>`__ in BioSTEAM, which covers most of the API, including unit operations. It also makes sure that the biorefineries are working properly. You can also just run each test separately for debuging. For example, the following runs the biorefinery tests:

.. code-block:: python
    
   >>> from biosteam.tests import test_biorefineries
   >>> test_biorefineries()
    

Changes made to a BioSTEAM unit operation requires it's specific unit test file to pass as well. For example:

.. code-block:: python

   >>> from biosteam.tests import test_binary_distillation
   >>> test_binary_distillation()

If no tests are available specfic to the unit operation, tests must be uploaded whereby the stream results and general simulation results are tested. Using `doctests <https://docs.python.org/3.6/library/doctest.html>`__ is the preferred method for running tests, but assertions are also accepted so long as all results of the unit operation is tested. You can use the source code in the "biosteam.tests" folder as a template for how this may be done.

.. note:: Several sections in biosteam are not fully documented. Any contributions towards rigorous testing is welcome!

Documentation
-------------

Concise and thorough documentation is required for any contribution. Make sure to:

* Use NumPy style docstrings.
* Document all functions and classes.
* Document short functions in one line if possible.
* Mention and reference any equations or methods used and make sure to include the chapter and page number if it is a book or a long document.
* Include a text file with the sphix autodoc in the "docs" folder.
* Preview the docs before making a pull request (open your cmd/terminal in the "docs" folder, run "make html", and open "docs/_build/html/index.html").
    
Best practices
--------------

Please refer to the following guides for best practices to make software designs more understandable, flexible, and maintainable:
    
* `PEP 8 style guide <https://www.python.org/dev/peps/pep-0008/>`__.
* `PEP 257 docstring guide <https://www.python.org/dev/peps/pep-0257/>`__.
* `Zen of Python philosophy <https://www.python.org/dev/peps/pep-0020/>`__.
* `SOLID programing principles <https://en.wikipedia.org/wiki/SOLID>`__.
