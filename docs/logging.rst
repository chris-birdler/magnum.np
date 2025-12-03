:tocdepth: 1

.. currentmodule:: magnumnp

#######
Logging
#######


Within *magnum.np* you are provided with three Loggers in order to log the simulation state. These are the :class:`.ScalarLogger`, the :class:`.FieldLogger` and the :class:`.Logger`. The :class:`.ScalarLogger` is useful to log scalar values such as the time and the averaged magnetization. The :class:`.FieldLogger` is used to log arbitrary scalar and vector fields alongside the simulation time. Lastly, :class:`.Logger` logs both scalar and field values by using both the :class:`.ScalarLogger` and the :class:`.FieldLogger`. Furthermore, *magnum.np* includes a resume function which enables the user to restart a simulation from the last logged state.

Logger
######
.. autoclass:: Logger
   :members:
   :special-members:
   :show-inheritance:
   :inherited-members:

The resume function is currently only available for the time *t* and the magnetization *m*. It uses existing log files and needs the field logger for the magnetization. Since this might not be given, resuming a state is not always possible and *False* otherwise. Therefore, you can check whether the state can be resumed by using :func:`.is_resumable()`. This function will return *True* in case resume is possible. Applying :func:`.resume()` will update the state object with the latest values of *t* and *m*. An example follows to show how this might be implemented. Here the magnetization is logged, then it is checked whether the state is resumable, after which the state is resumed and finally the state is logged.

.. code-block:: python

  logger = Logger('data', scalars = ['t', 'm'], fields = ['m'])
  logger.is_resumable()
  logger.resume(state)
  logger << state

ScalarLogger
############
.. autoclass:: ScalarLogger
   :members:
   :special-members:
   :show-inheritance:
   :inherited-members:

To explain how to log a custom scalar function we will look at further examples. In principle you can use :class:`.ScalarLogger` to log any function that accepts the state as a parameter and returns a scalar or a list/tuple. For example, you can log the state at a certain point. This would have to be done by defining a lambda function like so:

.. code-block:: python

  myfunc = lambda state: state.m[1,2,3,:]

This function would return the state at point [1,2,3,:]. To log this you can simply add *myfunc* to the list of variables that are being logged:

.. code-block:: python

  ScalarLogger(['m', 't', myfunc])

You might also want to log the average magnetization which you can do as follows:

.. code-block:: python

  myfunc = lambda state: state.m.avg()
  ScalarLogger(['m', 't', myfunc])

Or simply:

.. code-block:: python

    ScalarLogger(['m', 't', lambda state: state.m.avg()])


FieldLogger
###########
.. autoclass:: FieldLogger
   :members:
   :special-members:
   :show-inheritance:
   :inherited-members:

Just like the :class:`.ScalarLogger` the :class:`.FieldLogger` can also log custom functions. All you need is a function that takes the state as a parameter and returns a vectorfeld (nx,ny,nz,3). The name of the function is then simply given to the :class:`.FieldLogger` as a string. For example, to log the function *nsk* you would do the following:

.. Code-block:: python

  def nsk(state): # TODO: document and improve interface
      m = state.m.mean(axis=2)
      dxm = torch.stack(torch.gradient(m, spacing = state.mesh.dx[0], dim = 0), dim = -1).squeeze(-1)
      dym = torch.stack(torch.gradient(m, spacing = state.mesh.dx[1], dim = 1), dim = -1).squeeze(-1)
      return 1./(4.*pi) * (m * torch.linalg.cross(dxm, dym)).sum() * state.mesh.dx[0] * state.mesh.dx[1]

  FieldLogger("fields.pvd", ['m'], 'nsk', every  = 100)

Here *every = 100* was added, to log only every 100th step.
