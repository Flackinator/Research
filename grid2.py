import numpy as np
from matplotlib import pylab as plt

from physconst import AU, YR
from rotation import deg2rad
from color import rgb
from human import time2human

from multistar.generic import multi
from multistar.config import Config
from multistar.parallel import ParallelProcessor
from multistar.util import firsttrue

from multistar.interface import STATUS_OK, STATUS_COLLIDE, STATUS_ESCAPE

from multistar.grid.base import StudyBase, FateBase, TimeOutcome, SystemBase

class Fate(FateBase):
    FAIL = FateBase.FAIL
    STABLE = FateBase.STABLE
    UNKNOWN = FateBase.UNKNOWN
    COLLIDE = 10
    EARTHGONE = 20
    MOONGONE = 21

    labels = {
        FAIL: 'fail',
        STABLE: 'stable',
        COLLIDE: 'collision',
        EARTHGONE: 'earth gone',
        MOONGONE: 'moon gone',
        UNKNOWN: 'unknown',
        }

    keys = {v:k for k,v in labels.items()}

    colors = {
        FAIL : rgb('gray'),
        STABLE: rgb('k'),
        COLLIDE: rgb('r'),
        EARTHGONE: rgb('b'),
        MOONGONE: rgb('y'),
        UNKNOWN: rgb('g'),
        }

    assert np.all(np.array(list(colors.keys())) >= 0)

    colarr = np.zeros((np.max(list(colors.keys()))+1, 3))
    for i,v in colors.items():
        colarr[i,:] = v


class Outcome(TimeOutcome):
    fate = Fate


class Quad(SystemBase):

    vars = {
        'q' : ('q',  'mass ratio of binary', '$q={:5g}$'),
        'a' : ('an', 'semimajor axis of binary (AU)', '$a={:5g}$ (AU)'),
        'e' : ('en', 'eccentricity of binary', '$e={:5g}$'),
        'i' : ('i',  'inclination of binary (degrees)', '$i={:5g}$ (deg)'),
        'm' : ('pm', 'phase of the moon', r'$\phi_{{\mathrm{{moon}}}}={:5g}$ (deg)'),
        'b' : ('pb', 'phase of the binary', r'$\phi_{{\mathrm{{planet}}}}={:5g}$ (deg)'),
        }

    def __init__(self, toml='^/moon3.toml'):
        self.config = Config(toml)

    def __call__(self, en=0, an=0.1, i=0, q=1, pm=0, pb=0, dt=1*YR, cutoff=[0.01*AU, 0., 10*AU]):
        config = self.config.copy()
        if en is not None:
            config['binary orbit.en'] = en
        if an is not None:
            config['binary orbit.an_AU'] = an
        if i is not None:
            config['binary orbit.inclination_deg'] = i
        if q is not None:
            assert 0.1 <= q <= 1
            m1, m2 = np.array([1, q]) / (1 + q)
            config['Star A.M_Msun'] = m1
            config['Star B.M_Msun'] = m2

            # TODO - FIX: need to adjust radii properly
            config['Star A.S_Rsun'] = np.sqrt(m1)
            config['Star B.S_Rsun'] = np.sqrt(m2)
        if pm is not None:
            config['moon orbit.phase'] = pm
        if pb is not None:
            config['binary orbit.phase'] = pb
        if cutoff is not None:
            config.set('cutoff', cutoff)
        self.cutoff = cutoff

        # tx = np.minimum(dt, 10*YR)
        tx = dt
        dtd = tx
        return self.loop(config, dt, tx, dtd=dtd)

    def analyze(self, m, dt):

        # analize result

        if not hasattr(m, 't'):
            m.t = None
        if m.t is None:
            return Outcome(Fate.FAIL, 0)

        # test whether moon is further from earth than (solar) Hill radius (0.0098 AU)
        # ro = m.rom
        # if (ii := firsttrue(ro[0, :] > 0.01 * AU)) >= 0:
        #    return Outcome(Fate.MOONGONE, m.t[ii])
        if m.status == STATUS_ESCAPE and m.status_stars[0] == 1:
            return Outcome(Fate.MOONGONE, m.t[-1])

        # test whether moon and earth escaped jointly
        # since we already test whether moon escaped
        # if ((ii := np.argmax(ro[2, :] > self.cutoff)) > 0):
        #     return Outcome(Fate.EARTHGONE, m.t[ii] )
        # alternatively:
        if m.status == STATUS_ESCAPE:
            return Outcome(Fate.EARTHGONE, m.t[-1])

        if m.status == STATUS_COLLIDE:
            return Outcome(Fate.COLLIDE, m.t[-1])

        # unexpected case
        if m.status != STATUS_OK:
            return Outcome(Fate.UNKNOWN, m.t[-1] )

        return Outcome(Fate.STABLE, m.t[-1])


class Study(StudyBase):
    task = Quad
    fate = Fate
