import torch
import os
from magnumnp.common import logging
from xml.etree import ElementTree
from .io import write_vti 

__all__ = ["FieldLogger"]

class FieldLogger(object):
    def __init__(self, filename, fields, every = 1):
        """
        Logger class to fields

        *Arguments*
            filename (:class:`str`)
                The name of the log file
            fields ([:class:`str` | :class:`function`])
                The columns to be written to the log file
            every (:class:`int`)
                Write row to log file every nth call

        *Example*
            .. code-block:: python

                # provide key strings with are available in state
                logger = FieldLogger('data/m.pvd', ['m', demag.h])

                # Actually log fields
                state = State(mesh)
                logger << state
        """
        # create directory if not existent
        if not os.path.dirname(filename) == '' and \
             not os.path.exists(os.path.dirname(filename)):
            try:
                os.makedirs(os.path.dirname(filename))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        self._filename = filename
        self._file  = None
        self._every = every
        if isinstance(fields, str):
            fields = [fields]
        self._fields = fields
        self._i = 0

    def log(self, state):
        if self._file is None: # TODO: write pvd file
            self._file = open(self._filename+".pvd", 'w')
        self._i += 1
        if ((self._i - 1) % self._every > 0):
            return

        values = {}
        for field in self._fields:
            if isinstance(field, str):
                name = field
                value = getattr(state, field)
            elif hasattr(field, '__call__'):
                try:
                    name = field.__self__.__class__.__name__ + "." + field.__name__     
                except:
                    name = 'unnamed'
                value = field(state)
            elif isinstance(field, tuple) or isinstance(field, list):
                name = field[0]
                value = field[1](state)
            else:
                raise RuntimeError('Column type not supported.')
            values[name] = value

        write_vti(values, "%s_%04d.vti" % (self._filename, self._i // self._every), state = state)

    def __lshift__(self, state):
        self.log(state)

    def reset(self):
        self._i = 0

    def __del__(self):
        if self._file is not None:
            self._file.close()
