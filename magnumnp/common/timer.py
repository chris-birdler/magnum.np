#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import time
from collections import OrderedDict
from magnumnp.common import tabulate, logging
from functools import wraps

try:
    import resource
except:
    resource = None

__all__ = ["Timer", "timedmethod", "TimedOperator"]

class Timer(object):
    _timers    = OrderedDict()
    _start = None
    _current = ''
    _options = {
        'active': False,
        'skip': 0,
        'log_mem': False
    }

    def __init__(self, name):
        r"""
        Utility class that supports basic timing/profiling features. Each timer is
        initialized with a unique name. Each call to the body is timed and counted.
        Results can be printed with by calling :code:`print_report()`.

        *Example*
            .. code-block:: python

                Timer.enable()

                x = 0.0
                for i in range(10):
                    with Timer("Outer Function"):
                        for j in range(10):
                            with Timer("Inner Function"):
                                x += 1.0

                Timer.print_report()

            results in:

            .. code-block:: none

                ==================================================================
                TIMER REPORT
                ==================================================================
                Operation                     No of calls        Avg time [ms]        Total time [s]
                ----------------    -------------    ---------------    ----------------
                Outer Function                             10            0.0245333                 0.000245333
                    Inner Function                        100            0.000665188             6.65188e-05
                ----------------    -------------    ---------------    ----------------
                Total                                                                                                    0.00105381
                ==================================================================
        """
        if not Timer._options['active']: return

        self.t = 0.0
        self._fullname = Timer._current + '###' + name
        if self._fullname in self._timers:
            self._data= Timer._timers[self._fullname]
        else:
            self._data = {
                'name': name,
                'calls': 0,
                'total_time': 0.0,
                'memory': 0.0,
                'parent': Timer._current
            }
            Timer._timers[self._fullname] = self._data

    def __enter__(self):
        if not Timer._options['active']: return self

        Timer._current = self._fullname
        self._data['start'] = time.perf_counter()
        if Timer._options['log_mem']:
            self._data['mem_start'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return self

    def __exit__(self, e_typ, e_val, trcbak):
        if all((e_typ, e_val, trcbak)):
            raise e_typ from e_val

        if not Timer._options['active']: return self

        if Timer._options['skip'] <= self._data['calls']:
            self.t = time.perf_counter() - self._data['start']
            self._data['total_time'] += self.t
            if Timer._options['log_mem']:
                self.mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - self._data['mem_start']
                self._data['memory'] += self.mem
        self._data['calls'] += 1
        Timer._current = self._data['parent']
        return self

    @staticmethod
    def print_report():
        """
        Print the timing report.
        """
        def get_entries(parent = '', indent = ''):
            result = []

            for entry in filter(lambda x: x[1]['parent'] == parent, Timer._timers.items()):
                timed_calls = entry[1]['calls'] - Timer._options['skip']
                if timed_calls == 0:
                    t_avg = None
                    t_tot = None
                    if Timer._options['log_mem']:
                        mem = None
                else:
                    t_avg = (entry[1]['total_time'] / timed_calls) * 1000
                    t_tot = entry[1]['total_time']
                    if Timer._options['log_mem']:
                        mem = entry[1]['memory'] / 1024.

                if Timer._options['log_mem']:
                    result.append([indent + entry[1]['name'], entry[1]['calls'], t_avg, t_tot, mem])
                else:
                    result.append([indent + entry[1]['name'], entry[1]['calls'], t_avg, t_tot])
                result += get_entries(entry[0], indent + '    ')

            return result

        entries = get_entries()

        measured_time = sum([x[1]['total_time'] for x in Timer._timers.items() if x[1]['parent'] == ''])
        total_time = time.perf_counter() - Timer._start
        missing_time = total_time - measured_time

        if Timer._options['log_mem']:
            if Timer._start is not None:
                entries.append(['Total', None, None, total_time, None])
                entries.append(['Missing', None, None, missing_time, None])
            table = tabulate(entries, ["Operation", "No of calls", "Avg time [ms]", "Total time [s]", "Memory [MB]"])
        else:
            if Timer._start is not None:
                entries.append(['Total', None, None, total_time])
                entries.append(['Missing', None, None, missing_time])
            table = tabulate(entries, ["Operation", "No of calls", "Avg time [ms]", "Total time [s]"])

        # insert separator before total line
        lines = table.split('\n')
        lines.insert(len(lines) - 2, lines[1])
        table = "\n".join(lines)

        width = len(table.split("\n")[1])
        print("")
        print("=" * width)
        print("TIMER REPORT")
        print("=" * width)
        print(table)
        print("=" * width)
        print("")
        missing = missing_time / total_time
        if missing > 0.2:
            logging.warning("Too much time missing (%.0f%%). Add some Timers for more complete timing!" % (missing * 100.))

    @staticmethod
    def reset():
        """
        Reset all timers.
        """
        Timer._timers    = OrderedDict()
        Timer._current = ''
        Timer._options = {
            'active': False,
            'skip': 0,
            'log_mem': False
        }

    @staticmethod
    def configure(**kwargs):
        """
        Configure the timer.

        *Arguments*
            active :class:`bool`
                Flag the activates/deactived all timers.
            skip :class:`int`
                Skip the timing of the first x calls to every timer. The calls are
                still counted, but not considered in the average execution time.
            log_mem :class:`bool`
                Log accumulated memory consumption of timed calls
        """
        for (key, value) in kwargs.items():
            if key.lower() in Timer._options:
                Timer._options[key.lower()] = value
            else:
                raise ValueError("Option '%s' is not supported by Timer" % key)

        if Timer._options["log_mem"] == True and resource == None:
            raise RuntimeError("Import of 'resource' package failed. Try running Timer without 'log_mem = True'.")

    @staticmethod
    def disable():
        """
        Disable all timers.
        """
        Timer.configure(active = False)

    @staticmethod
    def enable(**kwargs):
        """
        Enable all timers. Takes all options that are accepted by :code:`configure`.
        """
        Timer._start = time.perf_counter()
        Timer.configure(active = True, **kwargs)


def timedmethod(method):
    @wraps(method)
    def timed(*args, **kwargs):
        with Timer("%s.%s" % (args[0].__class__.__name__, method.__name__)):
            return method(*args, **kwargs)
    return timed

def TimedOperator(obj, name="Unnamed Operator"):
    class wrapper(obj.__class__):
        def matvec(self, vec):
            with Timer(self._timer_name):
                return super(wrapper, self).matvec(vec)
    obj.__class__ = wrapper
    obj._timer_name = name
    return obj
