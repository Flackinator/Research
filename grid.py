
import numpy as np
from  matplotlib import pylab as plt

from physconst import AU, YR
from rotation import deg2rad
from color import rgb, ColorBlindRainbow
from human import time2human

from multistar.generic import multi
from multistar.config import Config
from multistar.parallel import ParallelProcessor
from multistar.util import firsttrue

from multistar.interface import STATUS_OK, STATUS_COLLIDE, STATUS_ESCAPE

from multistar.grid.base import StudyBase, FateBase, OutcomeBase, SystemBase



GOLDEN = 0.5 * (1 + np.sqrt(5))

class Fate(object):
    FAIL = 0
    STABLE = 1
    COLLISION = 10
    EARTHGONE = 20
    MOONGONE = 21

    labels = {
        FAIL: 'fail',
        STABLE: 'stable',
        COLLISION: 'collision',
        EARTHGONE: 'earth gone',
        MOONGONE: 'moon gone',
        }

    keys = {v:k for k,v in labels.items()}

    colors = {
        FAIL : rgb('gray'),
        STABLE: rgb('k'),
        COLLISION: rgb('r'),
        EARTHGONE: rgb('b'),
        MOONGONE: rgb('y'),
        }

    assert np.all(np.array(list(colors.keys())) >= 0)

    colarr = np.zeros((np.max(list(colors.keys()))+1, 3))
    for i,v in colors.items():
        colarr[i,:] = v

class Outcome(object):
    VERSION = 10000

    def __init__(self, outcome, time):
        if isinstance(outcome, str):
            outcome = Fate.keys[outcome]
        self.outcome = outcome
        self.time = time
        self.version = self.VERSION

    @property
    def stable(self):
        return self.outcome == Fate.STABLE

    @property
    def unstable(self):
        return self.outcome != Fate.STABLE

    def __str__(self):
        return f'{Fate.labels[self.outcome]}: {time2human(self.time)}'

    def __repr__(self):
        return f'{self.__class__.__name__}({str(self)})'

class Quad(object):

    def __init__(self, toml='binary_martin_base2.toml'):
        self.config = Config(toml)

    def __call__(self, en=0, an=0.1, i=0, q=1, pm=0, pb=0, dt=1*YR, cutoff=11*AU):
        config = self.config.copy()
        if en is not None:
            config['binary.2.en'] = en
        if an is not None:
            config['binary.2.an_AU'] = an
        if i is not None:
            config['binary.2.inclination_deg'] = i
        if q is not None:
            assert 0.1 <= q <= 1
            m1, m2 = np.array([1, q]) / (1 + q)
            config['star.3.M_Msun'] = m1
            config['star.4.M_Msun'] = m2

            # TODO - FIX: need to adjust radii properly - (DONE)
        if m1 <= 0.5:
            config['star.3.S_Rsun'] = m1**0.56
        else:
            config['star.3.S_Rsun'] = m1**0.79

        if m2 <= 0.5:
            config['star.4.S_Rsun'] = m2**0.56
        else:
            config['star.4.S_Rsun'] = m2**0.79

        if pm is not None:
            config['binary.1.phase'] = pm
        if pb is not None:
            config['binary.2.phase'] = pb
        if cutoff is not None:
            config.set('cutoff', cutoff)
        self.cutoff = cutoff


        m = multi(config)

        tx = np.minimum(dt, 1000*YR)
        tt = 0
        while True:
            m.rund(tx, dtd=0.1*YR)
            tt += tx
            outcome =  self.analyze(m, tt)
            if outcome.unstable:
                break
            if np.allclose(tt, dt):
                break
            tx = np.minimum(dt - tt, tx * GOLDEN)

        # for DEBUG only
        self.m = m
        self.config = config

        print(f' [{self.__class__.__name__}] {outcome!s}')

        return outcome

    def analyze(self, m, dt):

        # analize result

        if not hasattr(m, 't'):
            m.t = None
        if m.t is None:
            return Outcome(Fate.FAIL, 0)

        #if not np.allclose(m.t[-1], dt):
            # TODO - test what collided (use m.rn)
         #   return Outcome(Fate.COLLISION, m.t[-1])

        ro = m.ron
        # TODO - test whether moon is further from earth than (solar) hill radius
        #if ((ii := np.argmax(ro[0, :] > 0.01 * AU)) > 0):
           # return Outcome(Fate.MOONGONE, m.t[ii])
        # TODO - test whether moon and earth escaped jointly
        #if ((ii := np.argmax(ro[2, :] > 2 * AU)) > 0):
         #   return Outcome(Fate.EARTHGONE, m.t[ii] )

        if (ii := firsttrue(ro[0, :] > 0.01 * AU)) >= 0:
            return Outcome(Fate.MOONGONE, m.t[ii-1])
        if (ii := firsttrue(ro[2, :] > 10 * AU)) >= 0:
            return Outcome(Fate.EARTHGONE, m.t[ii-1])

        if not np.allclose(m.t[-1], dt):
            return Outcome(Fate.COLLISION, m.t[-1])

        return Outcome(Fate.STABLE, m.t[-1])


class Study(ParallelProcessor):
    VERSION = 10100

    vars = {
        'q' : ('q',  'mass ratio of binary', '$q={:5g}$'),
        'a' : ('an', 'semimajor axis of binary (AU)', '$a={:5g}$ (AU)'),
        'e' : ('en', 'eccentricity of binary', '$e={:5g}$'),
        'i' : ('i',  'inclination of binary (degrees)', '$i={:5g}$ (deg)'),
        'm' : ('pm', 'phase of the moon', r'$\phi_{{\mathrm{{moon}}}}={:5g}$ (deg)'),
        'b' : ('pb', 'phase of the binary', r'$\phi_{{\mathrm{{planet}}}}={:5g}$ (deg)'),
        }

    def legend(self, ax, rv = None):
        if rv is None:
            rv = list(Fate.colors.keys())
        for v in rv:
            ax.scatter(
                [None], [None],
                color=Fate.colors[v],
                label=Fate.labels[v],
                )
        leg = ax.legend(loc='best')
        leg.set_draggable(True)

    def __init__(self, **kwargs):
        kwargs.setdefault('task', Quad())
        super().__init__(**kwargs)

    def plot(self, vars=None, mode=None, data='fate'):
        fig, ax = plt.subplots()

        if vars is None:
            vars = ''
            for k,(v,l, _) in self.vars.items():
                vals = list()
                for r in self.results:
                    x =  getattr(r, v, None)
                    if x is not None:
                        if not x in vals:
                            vals.append(x)
                            if len(vals) > 1:
                                continue
                if len(vals) > 1:
                    vars += k
                    continue
        coords = list()
        labels = list()

        assert len(vars) == 2
        for k in vars:
            v, l, _ = self.vars[k]
            labels.append(l)
            coords.append(np.asarray([getattr(r, v) for r in self.results]))

        if data == 'time':
            rr = np.asarray([r.time for r in self.results.result()])
            c = np.log10(np.maximum(rr / YR, 1))
        elif data == 'fate':
            rr = np.asarray([r.outcome for r in self.results.result()])
            rv = np.unique(rr)
            c = Fate.colarr[rr]
        else:
            raise AttributeError(f'Unknown data "{data}".')

        if mode is None:
            dim = np.asarray([len(np.unique(cx)) for cx in coords])
            if np.product(dim) == c.shape[0]:
                mode = 'fill'

        if mode == 'fill':
            cb = list()
            ii = list()
            for cx in coords:
                cv = np.unique(cx)
                cc = np.ndarray(len(cv)+1)
                cc[1:-1] = 0.5 * (cv[1:] + cv[:-1])
                cc[[0,-1]] = 0.5 * (3 * cv[[0,-1]] - cv[[1,-2]])
                cb.append(cc)
                ii.append(np.searchsorted(cc, cx, side='right') - 1)

            z = np.ndarray((len(cb[0])-1, len(cb[1])-1) + c.shape[1:])
            z[ii[0], ii[1]] = c

            if data == 'time':
                cm = ax.pcolormesh(*cb, z.swapaxes(0,1), vmin=0, cmap=ColorBlindRainbow())
            else:
                ax.pcolorfast(*cb, z.swapaxes(0,1))
        else:
            ax.scatter(*coords, color=c)

        ax.set_xlabel(labels[0], fontsize=16)
        ax.set_ylabel(labels[1], fontsize=16)

        if data == 'time':
            ax==fig.colorbar(cm, label = 'log( time / years )')
        else:
            self.legend(ax, rv)

        fig.tight_layout()

    def __setstate__(self, state):
        super().__setstate__(state)
        if self.version < 10100:
            for r in self.results.results:
                r.result = Outcome(*r.result)
            self.version = 10100

    def __str__(self):
        return f'{self.__class__.__name__}(n={len(self.results)})'
