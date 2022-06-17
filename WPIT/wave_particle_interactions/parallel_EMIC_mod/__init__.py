import numpy as np
import environment_mod.const

from parallel_EMIC_mod.dzdt import dzdt
from parallel_EMIC_mod.dppardt import dppardt
from parallel_EMIC_mod.dpperdt import dpperdt
from parallel_EMIC_mod.detadt import detadt
from parallel_EMIC_mod.dlamdadt import dlamdadt
from parallel_EMIC_mod.dalphadt import dalphadt
from parallel_EMIC_mod.daeqdt import daeqdt
from parallel_EMIC_mod.dgammadt import dgammadt
from parallel_EMIC_mod.dEkdt import dEkdt
from parallel_EMIC_mod.nonlinear_S import nonlinear_S
from parallel_EMIC_mod.nonlinear_H import nonlinear_H
from parallel_EMIC_mod.nonlinear_theta import nonlinear_theta



