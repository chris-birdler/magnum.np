"""magnum.np main module"""

__version__ = '1.0.2'

try:
    import setproctitle
    setproctitle.setproctitle("magnumnp")

    from magnumnp.common import *
    from magnumnp.field_terms import *
    from magnumnp.solvers import *
    from magnumnp.loggers import *
    from magnumnp.utils import *

    import magnumnp.common.logging as logging
    logging.info_green("magnum.np %s" % __version__)

    import torch
    torch.set_default_dtype(torch.float64)

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
