import torch
import os
from collections.abc import Iterable
from magnumnp.loggers import ScalarLogger, FieldLogger
from magnumnp.common import logging

__all__ = ["Logger"]

class Logger(object):
    """
    Combined Scalar- and Field-Logger class

    *Arguments*
        directory (:class:`str`)
            The name of the log file
        scalars ([:class:`str` | :class:`function`])
            The columns to be written to the log file
        scalars_every (:class:`int`)
            Write scalar to log file every nth call
        fields ([:class:`str` | :class:`function`])
            The columns to be written to the log file
        every (:class:`int`)
            Write row to log file every nth call
        fields_every (:class:`int`)
            Write fields to log file every nth call

    *Example*
        .. code-block:: python

            # provide key strings with are available in state
            logger = Logger('data', ['m', demag.h], ['m'], fields_every = 100)

            # Actually log fields
            state = State(mesh)
            logger << state
    """
    def __init__(self, directory, scalars = [], fields = [], scalars_every = 1, fields_every = 1):
        self.loggers = {}
        if len(scalars) > 0:
            self.loggers["scalars"] = ScalarLogger(os.path.join(directory, "log.dat"), scalars, every = scalars_every)
        if len(fields) > 0:
            self.loggers["fields"] = FieldLogger(os.path.join(directory, "fields.pvd"), fields, every = fields_every)

    def log(self, state):
        for logger in self.loggers.values():
            logger << state

    def __lshift__(self, state):
        self.log(state)
