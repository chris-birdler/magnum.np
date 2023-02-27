"""magnum.np main module"""

__version__ = '1.0.8'

import magnumnp.common.logging as logging
import torch
#torch.set_default_dtype(torch.float64)

# monkey patch older torch version without compile
if not hasattr(torch, "compile"):
    def fake_compile(func):
        return func
    logging.warning("PyTorch version < 2.0 does not support 'torch.compile'! Switch to pytorch 2.0 in order get maximum performance!")
    torch.compile = fake_compile

try:
    import setproctitle
    setproctitle.setproctitle("magnumnp")

    from magnumnp.common import *
    from magnumnp.field_terms import *
    from magnumnp.solvers import *
    from magnumnp.loggers import *
    from magnumnp.utils import *

    logging.info_green("magnum.np %s" % __version__)

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
