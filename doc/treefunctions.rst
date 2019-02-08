Building expressions
====================

The :py:mod:`bamboo.treefunctions` module provides helper methods to construct expressions
and derived objects from the tree content (such that the user does not need to know about
the implementation details in the :py:mod:`bamboo.treeoperations` and :py:mod:`bamboo.treeproxies` modules).
In order to keep the code for such expressions compact, it is therefore recommended to import it with

.. code-block:: python

   from bamboo import treefunctions as op 

inside every analysis module.

List of functions
-----------------

.. automodule:: bamboo.treefunctions
   :members:
