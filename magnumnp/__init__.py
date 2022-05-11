"""magnum.np main module"""

VERSION = '0.0.1'

try:
    import setproctitle
    setproctitle.setproctitle("magnumnp")

    from torch import *

    from magnumnp.common import *
    from magnumnp.field_terms import *
    from magnumnp.solvers import *

    import magnumnp.common.logging as logging
    logging.info_green("magnum.np %s" % VERSION)

except Exception as e:
    import magnumnp.common.logging as logging
    logging.error(str(e).split("\n")[0])
    for line in str(e).split("\n")[1:]:
        logging.info(line)

    try:
        # do nothing if in IPython
        __IPYTHON__
        pass

    except NameError:
        # exit otherwise
        exit()
