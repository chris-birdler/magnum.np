"""magnum.np main module"""

__version__ = '1.0.9'

import magnumnp.common.logging as logging
import torch
torch.set_default_dtype(torch.float64)

# monkey patch torch versions < 2.0 or on Windows
try:
   def fake_compile(func):
       return func

   @torch.compile
   def dummy_function():
      return
   dummy_function()
except Exception as e:
    torch.compile = fake_compile
    logging.warning(e)


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
