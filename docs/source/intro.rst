.. raw:: html

    <script type="text/javascript"><!--
      function gen_mail_to_link(name, lhs,rhs)
      {
       document.write("<br>");
       document.write(name + " (<A HREF=\"mailto");
       document.write(":" + lhs + "@");
       document.write(rhs + "\">" + lhs + "@" + rhs + "<\/A>)"); } 
       // --> </SCRIPT>


About
=====

What is ALAMODE?
-----------------

**ALAMODE** is an open-source software package for analyzing lattice anharmonicity and lattice thermal conductivity of solids. Combined with an external DFT package such as VASP or Quantum ESPRESSO, ALAMODE extracts harmonic and anharmonic force constants in a straightforward way. The calculated anharmonic force constants can then be used to evaluate lattice thermal conductivity, phonon linewidths, and other anharmonic phonon properties from first principles.

Features
--------

The list of features, together with the force constants each of them needs, is given in :doc:`features`.

Links
-----

* Git repository : https://github.com/alamode-code/alamode (see :doc:`download`)
* Documentation  : https://alamode-code.github.io/alamode/ (this page)


Versioning
----------

Starting with 2.0.0, ALAMODE follows Semantic Versioning (MAJOR.MINOR.PATCH). Incompatible changes to input tags, file formats, or the Python API increase MAJOR; new features increase MINOR; bug fixes increase PATCH.


License
-------

.. |copy|   unicode:: U+000A9 

Copyright |copy| 2014-2026 Terumasa Tadano

This software is distributed under the MIT license.
See the LICENSE.txt file for license rights and limitations. 


How to Cite ALAMODE
-------------------

Please cite the appropriate articles from the list below when you publish results obtained with ALAMODE:

.. list-table:: ALAMODE papers
    :header-rows: 1

    * - Reference number
      - Features
    * - [Tadano2014]_
      - General citation paper for all features of ALAMODE
    * - [Tadano2015]_
      - Self-consistent phonon (SCP) method (``MODE = SCPH``)
    * - [Oba2019]_
      - Anharmonic free energies within SCP or improved SCP (ISC)
    * - [Masuki2022]_
      - SCP-based structural optimization or ``PERIODIC_IMAGE_CONV = 1`` (default) in the **alm** code
    * - [Masuki2023]_
      - QHA-based structural optimization

.. only:: html

   .. rubric:: References

.. [Tadano2014] T\. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter **26**\ , 225402 (2014) 
  [`Link <http://iopscience.iop.org/0953-8984/26/22/225402/>`__].

.. [Tadano2015] T\. Tadano and S. Tsuneyuki, Phys. Rev. B **92**\ , 054301 (2015). 
  [`Link <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.054301>`__]

.. [Oba2019] Y\. Oba, T. Tadano, R. Akashi, and S. Tsuneyuki, Phys. Rev. Materials **3**\, 033601 (2019). 
  [`Link <https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.3.033601>`__]

.. [Masuki2022] R\. Masuki, T. Nomoto, R. Arita, and T. Tadano, Phys. Rev. B **106**\, 224104 (2022).
  [`Link <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.224104>`__]

.. [Masuki2023] R\. Masuki, T. Nomoto, R. Arita, and T. Tadano, Phys. Rev. B **107**\, 134119 (2023).
  [`Link <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.134119>`__]


Issues & Bug Reports
--------------------

* If you find a bug or an issue related to ALAMODE, please report it at `GitHub issues <https://github.com/alamode-code/alamode/issues/>`_.
* Other questions and suggestions can be posted on the `GitHub discussions page <https://github.com/alamode-code/alamode/discussions>`_.

In either case, please search the existing posts and read the :ref:`FAQ <label_faq>` page of this document before asking.

Acknowledgment
--------------

This project has been partially supported by the following grants:

* JST PRESTO "Advancing first-principles prediction of physical properties through the introduction of finite-temperature effects"
* Grant-in-Aid for Scientific Research (C) (21K03424)
* Grant-in-Aid for Young Scientists (B) (16K17724) 
* Grant-in-Aid for Scientific Research on Innovative Areas 'Materials Design through Computics: Complex Correlation and Non-Equilibrium Dynamics'. (http://computics-material.jp)


Author & Contact
----------------

.. raw:: html

    <script>gen_mail_to_link('Terumasa TADANO', 'terumasa.tadano','gmail.com')</script>

| Research Center for Magnetic and Spintronic Materials (CMSM),
| National Institute for Materials Science (NIMS),
| Japan


Contributors
------------

* Ryota Masuki
* Tatsuro Nishimoto
* Yusuke Oba
* Yuto Tanaka
* Atsushi Togo
* Wenhao Zhang
