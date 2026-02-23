|pipeline status| |coverage report| |pypi| |conda| |version|

=====
Bilby
=====

A user-friendly Bayesian inference library.
Fulfilling all your Bayesian dreams.

Online material to help you get started:

-  `Installation instructions <https://bilby-dev.github.io/bilby/installation.html>`__
-  `Documentation <https://bilby-dev.github.io/bilby/>`__

If you need help, find an issue, or just have a question/suggestion you can

- Join our `Slack workspace <https://bilby-code.slack.com/>`__ via this `invite link <https://join.slack.com/t/bilby-code/shared_invite/zt-2s5a0jy1g-xB7uIy1fGGxW0CkBwwTbWQp>`__
- Ask questions (or search through other users questions and answers) on `StackOverflow <https://stackoverflow.com/questions/tagged/bilby>`__ using the bilby tag
- Submit issues directly through `the issue tracker <https://github.com/bilby-dev/bilby/issues>`__
- For chat.ligo.org users, join the `#bilby-help <https://chat.ligo.org/ligo/channels/bilby-help>`__ or `#bilby-devel <https://chat.ligo.org/ligo/channels/bilby-devel>`__ channels
- For LVK-confidential issues, please open `a confidential issue on bilby_pipe <https://git.ligo.org/lscsoft/bilby_pipe/-/issues/new>`__

We encourage you to contribute to the development of bilby. This is done via pull request.  For
help in creating a pull request, see `this page
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`__ or contact
us directly. For advice on contributing, see `the contributing guide <https://github.com/bilby-dev/bilby/blob/main/CONTRIBUTING.md>`__.


--------------
Citation guide
--------------

Please refer to the `Acknowledging/citing bilby guide <https://bilby-dev.github.io/bilby/citing-bilby.html>`__.

.. |pipeline status| image:: https://github.com/bilby-dev/bilby/actions/workflows/unit-tests.yml/badge.svg
   :target: https://github.com/bilby-dev/bilby/commits/master
.. |coverage report| image:: https://codecov.io/gh/bilby-dev/bilby/graph/badge.svg?token=BpfcD0hFSu
   :target: https://codecov.io/gh/bilby-dev/bilby
.. |pypi| image:: https://badge.fury.io/py/bilby.svg
   :target: https://pypi.org/project/bilby/
.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/bilby.svg
   :target: https://anaconda.org/conda-forge/bilby
.. |version| image:: https://img.shields.io/pypi/pyversions/bilby.svg
   :target: https://pypi.org/project/bilby/


-------------------------------------------
Instructions to use the tidalheating branch
-------------------------------------------

Local installation
~~~~~~~~~~~~~~~~~~

1. First, create and activate a conda environment (you can choose any name in place of "bilby-tidal"):

   ``conda create -n bilby-tidal python=3.10``

   ``conda activate bilby-tidal``

2. Install numpy<1.25, scipy<1.11 and contourpy<1.2 beforehand. At present more recent versions create conflicts. The order of installation of these packages is important.

   ``pip install "numpy<1.25"``

   ``pip install "scipy<1.11"``

   ``pip install "contourpy<1.2"``

3. Clone the repo (let's say in the home directory): 
   
   via https (works without ssh keys)

   ``git clone https://github.com/Samanwaya1301/bilby.git`` 
   
   or via SSH 

   ``git clone git@github.com:Samanwaya1301/bilby.git``

4. Enter and switch to the tidalheating branch:

   ``cd ~/bilby``

   ``git switch tidalheating``

5. Install dependencies:

   ``pip install -r requirements.txt``

   ``pip install -r gw_requirements.txt``

6. With an "editable install" of bilby, you won't have to install the package every time you make any change in the source repo. The package will be loaded directly from the source. For this, run:

   ``pip install -e .``

   With a normal installation instead, the package has to be reinstalled every time a source code change is expected to show up in the outputs. For this, run (and run this after every change): 

   ``pip install .``

   Each of these options has its pros and cons, and the user should make the choice considering their use case.

7. Test the installation:

   ``cd ~/bilby/examples/gw_examples/injection_examples``

   ``python htf2_check.py``


Installation in LDG clusters
~~~~~~~~~~~~~~~~~~~~~~~~

For a cluster installation, *editable install is not recommended* as running jobs may take several days and the results may be inconsistent if the source code is changed during that period.

If the cluster is equipped with IGWN conda distributions and you choose to use that, you may create a conda environment with 

``conda create -n bilby-tidal --clone igwn-py310``

and install the repo by cloning and activating the tidalheating branch. This will override the already existing Bilby installation in the IGWN clone.

If you don't clone IGWN environment and follow the above steps for a minimal installation in the cluster, you may want to also install `bilby_pipe <https://lscsoft.docs.ligo.org/bilby_pipe/0.3.12/index.html>`_ for job submissions: 

``conda install -c conda-forge bilby_pipe``

