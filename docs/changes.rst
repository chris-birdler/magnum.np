.. module:: magnumnp

:tocdepth: 1

########################
Changes in magnum.np 2.0
########################

During the first major refactoring of the `magnum.np` code the following problems were identified:

* ``torch.compile`` does not work with ``DecoratedTensors`` (in fact it works, but it gives bad timing since it cannot optimize the code correctly)
* rarely used features like non-equidistant meshes influenced the main-stream code 
* ``torch.Tensor``'s directly created by user may end up on wrong `device` or with wrong `dtype`

The following design changes have been performed to address these issues:

************************
Remove Decorated Tensors
************************

``DecoratedTensors`` are completely removed in order to allow directly using ``torch.compile``. Previous `magnum.np` versions used wrapper functions which used only the raw ``torch.Tensor``. Only those low-level functions could be compiled efficiently. However, this workaround complicates the code and prevents the compilations of larger code entities (like e.g. the full LLG integrator including all effective field terms).

In order to remove the ``DecoratedTensor`` the following features have been re-implemented:

* ``DecoratedTensor.normalize()`` method has been moved to ``normalize(tensor)`` function. 
* ``DecoratedTensor.avg()`` method has been moved to ``State.avg(tensor)``. This is necessary because averaging on non-equidistant meshes requires the knowledge of the cell volumes, in order to use correct weights. 
* ``DecoratedTensor.__call__()`` was originally used to generalize ``lambda`` functions and constant tensors. Internaly all material parameters could be used as if they were ``lambda`` functions. In case of the ``DecoratedTensor`` its call function simply returned the tensor itself. In `magnum.np 2.0` this mechanism has been re-implemented in the Material setter. When a constant material is set, a dummy lambda function is created which returns the constant tensor. Similar checks have to be done for ``State.j`` and ``ExternalField.h``.

**********************
Non-Equidistant Meshes
**********************
In `magnum.np 2.0` the ``Mesh`` object contains a ``is_equidistant`` flag which tells whether it is equidistant or not. There is an extended tensor ``dx`` which allow a generalization of exchange-like field terms as well as an ``dx_tuple`` containing the original grid-spacing as a tuple (this is used for code which has not been generalized to non-equidistant mesh). 

The most important change is the handling of slices. In previous versions each tensor carried a copy of the ``cell_volumes`` tensor with itself in order to allow calculation of a proper average if needed. This required additional memory and also led to problems which slices of a tensor were used (in this case also the ``cell_volumes`` had to be sliced with a modified slice). 

In `magnum.np 2.0` the average method has been moved to the ``State``. In case of a equidistant mesh is calculate a simple algebraic mean. In case of non-equidistant meshes it checks whether a full-tensor or a slice is provided. In case of a slice the user needs to provice the sliced ``cell_volumes`` as additional parameter. Otherwise an Exception will be raised. In this way the traditional logging of full tensors directly works and the user gets a proper error message for more complicated scenarios.

************************
Setting Device and dtype
************************
In `magnum.np 2.0` the `device` and `dtype` are set by changing the PyTorch defaults (an optional State parameter is no longer supported). This allows to remove the wrapper functions for Pytorch generator functions like ``arange`` or ``linspace`` inside of the ``State``. The user can directly use pytorch functions and the tensors will be created on the correct `device` with the correct `dtype`.


****************
Lambda Functions
****************
All lambda functions are now depending on the ``State`` and no longer only on the time. This generalizes the material dependency and allows e.g. also temperature or current dependent parameters. Unfortunately, it requires some slight modification of the user code e.g.:

.. code-block:: python

    # old code
    h_ext = ExternalField(lambda t: [0, 0, t])

    # new code
    h_ext = ExternalField(lambda state: [0, 0, state.t])


