# SPDX-License-Identifier: MIT
"""Independent two-dimensional metrics for the WJ normal hill."""

import math
import numpy as np
from scipy.optimize import brentq


class Hill:
    def __init__(self, geometry):
        self.left, self.right = geometry["left"], geometry["right"]
        self.H, self.d = geometry["H"], geometry["d0"]
        self.top, self.height = geometry["tunnel_top"], geometry["hill_height"]
        self.L = geometry["hill_half_length"]
        self.width = self.top - self.d - self.H

    def graph(self, a, foot_side=0):
        a = np.asarray(a)
        w = np.pi / self.L
        at_foot = np.abs(np.abs(a) - self.L) <= 32 * np.finfo(float).eps * max(1., self.L)
        interior = (np.abs(a) < self.L) & ~at_foot
        inside_limit = at_foot & (a * foot_side < 0)
        f = np.where(interior, self.height * np.cos(w * a / 2)**2, 0.)
        p = np.where(interior, -self.height * w / 2 * np.sin(w * a), 0.)
        pp = np.where(interior | inside_limit, -self.height * w**2 / 2 * np.cos(w * a), 0.)
        return f, p, pp

    def metrics(self, e, a, r, outer=None):
        xl, xr, rl, _ = e["bounds"]
        foot_side = np.where(np.abs(a - xl) <= 2e-14, 1,
                             np.where(np.abs(a - xr) <= 2e-14, -1, 0))
        f, p, pp = self.graph(a, foot_side)
        ell = np.hypot(1., p)
        tangent = np.stack((1 / ell, p / ell), axis=-1)
        normal = np.stack((-p / ell, 1 / ell), axis=-1)
        kappa = -pp / ell**3
        kh = kappa / (1 + (self.d + self.H) * kappa)
        if outer is None:
            outer = rl >= self.H - 1e-13
        if outer:
            d = self.d + self.H
            xa, za = a + d * normal[:, 0], f + d * normal[:, 1]
            eta = (r - self.H) / self.width
            position = np.column_stack(((1 - eta)*xa + eta*a, (1 - eta)*za + eta*self.top))
            along = ((1 - eta)*(ell*(1 + d*kappa)))[:, None]*tangent
            along[:, 0] += eta
            radial = np.column_stack((a - xa, self.top - za)) / self.width
        else:
            d = self.d + r
            position = np.column_stack((a, f)) + d[:, None]*normal
            along = (ell*(1 + d*kappa))[:, None]*tangent
            radial = normal
        det = along[:, 0]*radial[:, 1] - along[:, 1]*radial[:, 0]
        if np.any(det <= 0):
            raise ValueError("independent chart Jacobian failed")
        da = np.column_stack((radial[:, 1], -radial[:, 0])) / det[:, None]
        dr = np.column_stack((-along[:, 1], along[:, 0])) / det[:, None]
        return dict(position=position, det=det, da=da, dr=dr, kh=kh, t=tangent, n=normal,
                    effective_area=ell*(1 + self.d*kappa),
                    matching_area=ell*(1 + (self.d + self.H)*kappa),
                    curvature=kappa/(1 + self.d*kappa))

    def face(self, e, side, a, r):
        m = self.metrics(e, a, r)
        area = m["det"][:, None]*m["da" if side < 2 else "dr"]*(1 if side % 2 else -1)
        length = np.linalg.norm(area, axis=1)
        return length, area/length[:, None]

    def certificate(self):
        bound = self.height * (math.pi/self.L)**2 / 2
        reach = self.d + self.H
        slope = self.height * math.pi/(2*self.L)
        minimum_tube_stretch = 1 - reach*bound
        roof_gap = self.top - self.height - reach*math.hypot(1., slope)
        if minimum_tube_stretch <= 0 or roof_gap <= 0:
            raise ValueError("global chart certificate failed")
        # For a C1,1 graph this bound makes each point strictly nearer
        # its own physical hill foot than any other foot. For the ruled
        # outer mesh x_a>0 and z rises at fixed x by det/x_a.
        return dict(physical_graph_curvature_bound=bound,
                    minimum_tube_stretch=minimum_tube_stretch,
                    minimum_effective_J=minimum_tube_stretch/(1 + self.d*bound),
                    minimum_outer_determinant=minimum_tube_stretch*roof_gap/self.width)

    def from_experiment(self, x_mm, above_hill_mm):
        """Map a vertical laboratory sample at fixed x into the chart."""
        x = x_mm/40
        tolerance = 64*np.finfo(float).eps*max(1.,abs(self.left),abs(self.right),abs(self.top))
        if not self.left-tolerance <= x <= self.right+tolerance:
            raise ValueError("experiment sample is outside the tunnel mesh")
        # Unit conversion can place a boundary sample a few ulps past a face.
        x = min(self.right,max(self.left,x))
        z = float(self.graph(x)[0]) + above_hill_mm/40
        self.certificate()
        if not math.isfinite(z) or z > self.top+tolerance:
            raise ValueError("experiment sample is outside the tunnel mesh")
        z = min(z,self.top)
        reach = self.d + self.H
        limits = [x - reach, x + reach]
        def matching(a):
            f, p, _ = self.graph(a)
            ell = math.hypot(1., p)
            return a - reach*p/ell, f + reach/ell
        foot = brentq(lambda a: matching(a)[0] - x, *limits)
        if z <= matching(foot)[1]:
            def orthogonality(a):
                f, p, _ = self.graph(a)
                return a - x + (f - z)*p
            a = brentq(orthogonality, *limits)
            f, p, _ = self.graph(a)
            r = ((z - f) - (x - a)*p)/math.hypot(1., p) - self.d
            if r < -tolerance:
                raise ValueError("experiment sample is below the displaced wall")
            return a, max(0.,r)
        def projected_x(a):
            xh, zh = matching(a)
            eta = (z - zh)/(self.top - zh)
            return (1 - eta)*xh + eta*a
        a = brentq(lambda a: projected_x(a) - x, *limits)
        zh = matching(a)[1]
        return a, self.H + self.width*(z - zh)/(self.top - zh)
