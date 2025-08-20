"""magnum.np main module"""

__version__ = '2.1.0'

import magnumnp.common.logging as logging
from magnumnp.common.utils import get_gpu_with_least_memory
import torch
import os
import sys

# set default device and dtype
torch.set_default_dtype(torch.float64)
device_id = os.environ.get('CUDA_DEVICE')
if device_id == None:
    device_id = get_gpu_with_least_memory()
device = torch.device(f"cuda:{device_id}" if int(device_id) >= 0 else "cpu")
torch.cuda.set_device(int(device_id)) # prevents nvidia-smi from showing process on GPU0 due to torch.compile
torch.set_default_device(device)
torch.manual_seed(2147483647) # fix seed

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
    setproctitle.setproctitle("magnumnp %s" % sys.argv[0])

    from magnumnp.common import *
    from magnumnp.field_terms import *
    from magnumnp.linear_elasticity import *
    from magnumnp.solvers import *
    from magnumnp.loggers import *
    from magnumnp.utils import *

    logging.info_green("magnum.np %s (%s)" % (__version__, (" ".join(sys.argv))))

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
