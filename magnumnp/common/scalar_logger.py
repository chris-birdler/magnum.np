import torch
from collections.abc import Iterable
from functools import reduce
import os
from magnumnp.common import logging, DecoratedTensor 

__all__ = ["ScalarLogger"]

class ScalarLogger(object):
    def __init__(self, filename, columns, every = 1):
        """
        Simple logger class to log scalar values into a tab separated file.

        *Arguments*
            filename (:class:`str`)
                The name of the log file
            columns ([:class:`str` | :class:`function`])
                The columns to be written to the log file
            every (:class:`int`)
                Write row to log file every nth call

        *Example*
            .. code-block:: python

                # provide key strings with are available in state
                logger = ScalarLogger('log.dat', ['t','m'])

                # provide func(state) or tuple (name, func(state))
                logger = ScalarLogger('log.dat', [('t[ns]', lambda state: state.t*1e9)])

                # provide predifined functions
                logger = ScalarLogger('log.dat', [demag.h, demag.E])

                # Actually log a row
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
        self._file     = None
        self._every    = every
        self._i        = 0
        self._columns  = columns

    def add_column(self, column):
        if self._file is not None:
            raise RuntimeError("You cannot add columns after first log row has been written.")
        self._columns.append(column)

    def log(self, state):
        self._i += 1
        if ((self._i - 1) % self._every > 0):
            return

        values = []

        for column in self._columns:
            if isinstance(column, str):
                name = column
                raw_value = getattr(state, column)
            elif hasattr(column, '__call__'):
                try:
                    name = column.__self__.__class__.__name__ + "." + column.__name__     
                except:
                    name = 'unnamed'
                raw_value = column(state)
            elif isinstance(column, tuple) or isinstance(column, list):
                name = column[0]
                raw_value = column[1](state)
            else:
                raise RuntimeError('Column type not supported.')

            if isinstance(raw_value, DecoratedTensor):
                value = raw_value.avg().tolist()
            elif isinstance(raw_value, torch.Tensor):
                value = raw_value.tolist()
            else:
                value = raw_value
            values.append((name, value))
        if self._file is None:
            self._file = open(self._filename, 'w')
            self._write_header(values)

        self._write_row(values)

    def __lshift__(self, state):
        self.log(state) 

    def _write_header(self, columns):
        headings = []

        for column in columns:
            if isinstance(column[1], Iterable):
                if (len(column[1]) == 3):
                    for i in ('x', 'y', 'z'):
                        headings.append(column[0] + '_' + i)
                else:
                    for i in range(len(column[1])):
                        headings.append(column[0] + '_' + str(i))
            else:
                headings.append(column[0])

        format_str = "#" + "    ".join(["%-22s"] * len(headings)) + "\n"
        self._file.write(format_str % tuple(headings))
        self._file.flush()
        
    def _write_row(self, columns):
        flat_values = reduce(lambda x,y: x+y,
            map(lambda x: tuple(x[1]) if isinstance(x[1], Iterable) else (x[1],), columns))
        format_str = "    ".join(["%+.15e"] * len(flat_values)) + "\n"
        self._file.write(format_str % flat_values)
        self._file.flush()

    def __del__(self):
        if self._file is not None:
            self._file.close()
