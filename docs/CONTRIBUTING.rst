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

<<<<<<< HEAD
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
=======
First pip install `pytest <https://docs.pytest.org/en/stable/>`__ and run the
following in your local biosteam directory:

.. code-block:: bash
    
   $ pytest --doctest-modules
    =================================== test session starts ===================================
    platform win32 -- Python 3.7.6, pytest-5.3.5, py-1.8.1, pluggy-0.13.1
    rootdir: C:\Users\...\biosteam, inifile: pytest.ini
    plugins: hypothesis-5.5.4, arraydiff-0.3, astropy-header-0.1.2, doctestplus-0.5.0, 
    openfiles-0.4.0, remotedata-0.3.2
    collected 30 items
    
    biosteam\_heat_utility.py ..                                                         [  6%]
    biosteam\_power_utility.py .                                                         [ 10%]
    biosteam\examples\ethanol_subsystem_example.py .                                     [ 13%]
    biosteam\process_tools\unit_group.py .                                               [ 16%]
    biosteam\units\_balance.py .                                                         [ 20%]
    biosteam\units\_binary_distillation.py .                                             [ 23%]
    biosteam\units\_duplicator.py .                                                      [ 26%]
    biosteam\units\_fermentation.py .                                                    [ 30%]
    biosteam\units\_flash.py .                                                           [ 33%]
    biosteam\units\_hx.py ..                                                             [ 40%]
    biosteam\units\_junction.py .                                                        [ 43%]
    biosteam\units\_liquids_centrifuge.py .                                              [ 46%]
    biosteam\units\_lle_unit.py .                                                        [ 50%]
    biosteam\units\_mixer.py .                                                           [ 53%]
    biosteam\units\_molecular_sieve.py .                                                 [ 56%]
    biosteam\units\_process_specification.py .                                           [ 60%]
    biosteam\units\_pump.py .                                                            [ 63%]
    biosteam\units\_shortcut_column.py .                                                 [ 66%]
    biosteam\units\_splitter.py .                                                        [ 70%]
    biosteam\units\_tank.py .                                                            [ 73%]
    biosteam\units\design_tools\heat_transfer.py .                                       [ 76%]
    biosteam\units\design_tools\tank_design.py ..                                        [ 83%]
    biosteam\utils\functors.py .                                                         [ 86%]
    biosteam\utils\piping.py ..                                                          [ 93%]
    tests\test_example_sugarcane_subsystem.py ..                                         [100%]
    
    =================================== 30 passed in 5.84s ====================================
    
This runs all the `doctests <https://docs.python.org/3.6/library/doctest.html>`__
in BioSTEAM, which covers most of the API, including unit operations. If any test 
is marked with a letter F, that test has failed. Pytest will point you to the
location of the error, the values that were expected, and the values that were 
generated.

Changes made to a BioSTEAM unit operation requires it's specific doctests to pass
before uploading. If no tests are available specfic to the unit operation, tests 
must be uploaded whereby the stream results and general simulation results are 
tested. Using `doctests <https://docs.python.org/3.6/library/doctest.html>`__ is
the preferred method for running tests, but assertions in a test function is also
accepted so long as all results of the unit operation is tested. You can use the
source code in the "biosteam.tests" folder as a template for how this may be done.

.. note:: 

    Several sections in biosteam are not fully documented. Any contributions
    towards rigorous testing is welcome!

The `biorefineries <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park>`__ 
package can also be tested by running the following in your local biorefineries repository:

.. code-block:: bash

    $ pytest
    =================================== test session starts ===================================
    platform win32 -- Python 3.7.6, pytest-5.3.5, py-1.8.1, pluggy-0.13.1
    rootdir: C:\Users\...\Bioindustrial-Park\BioSTEAM 2.x.x
    plugins: hypothesis-5.5.4, arraydiff-0.3, astropy-header-0.1.2, doctestplus-0.5.0, 
    openfiles-0.4.0, remotedata-0.3.2
    collected 2 items
    
    tests\test_biorefineries.py ..                                                       [100%]
    
    =================================== 2 passed in 4.62s =====================================

The `thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`__ 
package is not yet ready for `pytest`. Instead, run the following lines in python to
test thermosteam:

.. code-block:: python

    >>> from thermosteam.tests import test_thermosteam  # From local repository
    >>> test_thermosteam() # If nothing happens, all tests have passed
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1

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
