import numpy as np

from physconst import AU, YR
from rotation import deg2rad

from multistar.config import Config
from multistar.util import firsttrue

from multistar.interface import STATUS_OK, STATUS_COLLIDE, STATUS_ESCAPE

from multistar.grid.base import StudyBase, FateBase, Outcome, SystemBase

class Fate(FateBase):
    FAIL = FateBase.FAIL
    STABLE = FateBase.STABLE
    COLLISION = 10
    EARTHGONE = 20
    MOONGONE = 21
    UNKNOWN = 99

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

class Quad(SystemBase):

    vars = {
        'q' : ('q',  'mass ratio of binary', '$q={:5g}$'),
        'a' : ('an', 'semimajor axis of binary (AU)', '$a={:5g}$ (mAU)'),
        'e' : ('en', 'eccentricity of binary', '$e={:5g}$'),
        'i' : ('i',  'inclination of binary (degrees)', '$i={:5g}$ (deg)'),
        'm' : ('pm', 'phase of the moon', r'$\phi_{{\mathrm{{moon}}}}={:5g}$ (deg)'),
        'b' : ('pb', 'phase of the binary', r'$\phi_{{\mathrm{{planet}}}}={:5g}$ (deg)'),
        }

    def __init__(self, toml='^/moon3.toml'):
        self.config = Config(toml)

    def __call__(self, en=0, an=0.1, i=0, q=1, pm=0, pb=0, dt=1*YR, cutoff=10*AU):
        config = self.config.copy()
        if en is not None:
            config['binary.2.en'] = en
        if an is not None:
            config['binary.2.an_AU'] = an
        if i is not None:
            euler_deg = config['binary.2.euler_deg']
            euler_deg[1] = i
            config['binary.2.euler_deg'] = euler_deg
        if q is not None:
            assert 0.1 <= q <= 1
            m1, m2 = np.array([1, q]) / (1 + q)
            config['star.3.M_Msun'] = m1
            config['star.4.M_Msun'] = m2

            # TODO - FIX: need to adjust radii properly
            config['star.3.S_Rsun'] = np.sqrt(m1)
            config['star.4.S_Rsun'] = np.sqrt(m2)
        if pm is not None:
            config['binary.1.phase'] = pm
        if pb is not None:
            config['binary.2.phase'] = pb
        if cutoff is not None:
            config.set('cutoff', cutoff)
        if cutoff is None:
            cutoff = 2 * AU
        self.cutoff = cutoff

        tx = np.minimum(dt, 10*YR)
        return self.loop(config, dt, tx)

    def analyze(self, m, dt):

        # analize result

        if not hasattr(m, 't'):
            m.t = None
        if m.t is None:
            return Outcome(Fate.FAIL, 0)

        ro = m.ron
        # test whether moon is further from earth than (solar) Hill radius (0.0098 AU)
        if (ii := firsttrue(ro[0, :] > 0.01 * AU)) >= 0:
            return Outcome(Fate.MOONGONE, m.t[ii])

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


class Study(ParallelProcessor):
    task = Quad
    fate = Fate
