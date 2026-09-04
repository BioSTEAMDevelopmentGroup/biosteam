Security
========
BioSTEAM is generally safe to use for process simulation and engineering 
analysis. As with any Python software, however, users should exercise 
appropriate care when working with untrusted inputs.

Security issues can be reported privately to the BioSTEAM maintainers through 
`GitHub <https://github.com/BioSTEAMDevelopmentGroup/biosteam/security/advisories?state=Triage>`_. 
Please do not disclose a suspected vulnerability through a public issue 
until the maintainers have had an opportunity to assess and address it.

When reporting a security issue, please provide enough information to 
reproduce and assess the problem, including the affected BioSTEAM version, 
relevant code or input data, and a description of the potential impact.

Supported versions
------------------
Security updates are provided for BioSTEAM 2.54.0 and later. 
Versions older than 2.54.0 are not supported and do not receive security updates.
Users of unsupported versions are encouraged to upgrade to a supported release.

Safety considerations
---------------------
BioSTEAM is generally safe to use with native models and user data. 
However, some Python functionality can execute code or evaluate expressions, 
and care should be taken when processing data or workflows from untrusted sources.
If BioSTEAM is used to process untrusted data or execute workflows provided 
by untrusted users, appropriate isolation or sandboxing is recommended. 

Reporting
---------
Please report security issues privately using 
`GitHub's Security Advisories feature <https://github.com/BioSTEAMDevelopmentGroup/biosteam/security/advisories?state=Triage>`_ 
rather than through the public issue tracker.
Security reports are reviewed by the BioSTEAM maintainers and, when 
appropriate, a fix will be developed and released through a supported 
BioSTEAM version.