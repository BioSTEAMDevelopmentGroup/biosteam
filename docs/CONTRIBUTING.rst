Contributing to BioSTEAM
========================

Monthly coordination calls
--------------------------

We welcome anyone with an interest in discussing contributions to BioSTEAM to join our
monthly meeting at 8-9am CDT, every 3rd Friday. Please email biosteamdevelopmentgroup@gmail.com 
for a link to join. Here is the agenda for this month's coordination call:

.. toctree::
   :maxdepth: 2
   
   coordination/November-19-2021

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

First pip install `pytest <https://docs.pytest.org/en/stable/>`__ and run the
following in your local biosteam directory:

.. code-block:: bash
    
   $ pytest
    
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
accepted so long as all results of the unit operation is tested. 

.. note:: 

    Several sections in biosteam are not fully documented. Any contributions
    towards rigorous testing is welcome!

The `biorefineries <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park>`__ 
package can also be tested by running pytest in your local biorefineries repository:

If you are uploading a new biorefinery to the `biorefineries` package, make sure
you include tests for the following results:

* One TEA feasibility parameter (e.g. MPSP, MFSP, or IRR).
* Sales
* Material cost
* Installed equipment cost
* Utility cost
* Heating duty
* Coling duty
* Electricity consumption
* Electricity production 

Also make sure to add a README.rst file to explain basic functionality and 
design. If a paper is already published, you can add a link to it here too.
Examples in the readme may be tested using `doctest's testfile <https://docs.python.org/3/library/doctest.html>`__
method. In the biorefineries repository you can find `example tests <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/tests>`__.

The `thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`__ 
package also uses `pytest` for testing (just run pytest on your local 
thermosteam repository).

Documentation
-------------

Concise and thorough documentation is required for any contribution. Make sure to:

* Use NumPy style docstrings.
* Document all functions and classes.
* Document short functions in one line if possible.
* Mention and reference any equations or methods used and make sure to include the chapter and page number if it is a book or a long document.
* Preview the docs before making a pull request (open your cmd/terminal in the "docs" folder, run "make html", and open "docs/_build/html/index.html").

Authorship
----------

Authorship must be acknowledged for anyone contributing code, significant  
expertise, and/or other impactful efforts. Additionally, authorship should be 
included at the module-level, with a short description of the general contribution. 

If any code or implementation was copied from a third party, it should be rightfully
noted in the module-level documentation along with the corresponding copyright.

Any third-party code copied to the BioSTEAM software must be strictly open-source 
(not copy-left nor open-access). Additionally, if the license is different, 
the module should add the third-party license as an option (dual licensing is OK).


Best practices
--------------

Please refer to the following guides for best practices to make software designs more understandable, flexible, and maintainable:
    
* `PEP 8 style guide <https://www.python.org/dev/peps/pep-0008/>`__.
* `PEP 257 docstring guide <https://www.python.org/dev/peps/pep-0257/>`__.
* `Zen of Python philosophy <https://www.python.org/dev/peps/pep-0020/>`__.
* `SOLID programing principles <https://en.wikipedia.org/wiki/SOLID>`__.
