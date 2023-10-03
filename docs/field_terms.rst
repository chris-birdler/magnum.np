.. module:: magnumnp

.. _Sphinx Fields:

###########
Field Terms
###########

The effective field can be comprised of many different fields such as the quantum mechanical exchange interaction or an anisotropy field due to spin orbit coupling. :class:`FieldTerm` is a class representing the possible field terms. magnum.np comes with a number of predefined field term classes that cover a variety of use cases. In most applications a list of :class:`FieldTerm` objects is usually passed to the solver in order to solve micromagnetic problems under the influence of different effects. Furthermore, each :class:`FieldTerm` object has methods to compute the corresponding effective field and energy for a given State object.

Here is an examples on how to use the field term classes:

.. code-block:: python

  # initialize state with uniform magnetization configuration
  state = State(mesh, m = Constant((1.0, 0.0, 0.0)))

  # initialize demagnetization field object
  demag = DemagField()

  # compute field and save to file
  File("h.pvd") << demag.h(state)

  # compute energy and print on screen
  print demag.E(state)

Both, the field method :code:`h` and the energy method :code:`E` can also be used as virtual attribues:

.. code-block:: python

  state.h_demag = demag.h()
  state.E_demag = demag.E()

The computation of the field is carried out in the default region which usually is the ‘magnetic’ region. However, if the user wishes to restrict the computation of the field to a certain region, this default can be overriden in the initialization, like so:

.. code-block:: python

  exchange = ExchangeField(region = 'ferromagnetic')

where the region called 'ferromagnetic' would have to be defined by the user beforehand. How this can be done is explained under :ref:`Sphinx Domains`.

Some field term classes require certain material parameters to be defined. The following documentatin offers more detailed information about each class.

LinearFieldTerm
***************
.. autoclass:: LinearFieldTerm
   :show-inheritance:

UniaxialAnisotropyField
***********************
.. autoclass:: UniaxialAnisotropyField
   :show-inheritance:

CubicAnisotropyField
********************
.. autoclass:: CubicAnisotropyField
   :show-inheritance:

DemagField
**********
.. autoclass:: DemagField
   :show-inheritance:

DemagFieldPBC
*************
.. autoclass:: DemagFieldPBC
   :show-inheritance:

DemagFieldNonEquidistant
*************
.. autoclass:: DemagFieldNonEquidistant
   :show-inheritance:

ExchangeField
*************
.. autoclass:: ExchangeField
   :show-inheritance:

ExternalField
*************
.. autoclass:: ExternalField
   :show-inheritance:

DMIField
********
.. autoclass:: DMIField
   :show-inheritance:

BulkDMIField
************
.. autoclass:: BulkDMIField
   :show-inheritance:

InterfaceDMIField
*****************
.. autoclass:: InterfaceDMIField
   :show-inheritance:

D2dDMIField
***********
.. autoclass:: D2dDMIField
   :show-inheritance:

RKKYField
*********
.. autoclass:: RKKYField

SpinOrbitTorque
***************
.. autoclass:: SpinOrbitTorque

SpinTorqueZhangLi
*****************
.. autoclass:: SpinTorqueZhangLi

OerstedField
************
.. autoclass:: OerstedField
