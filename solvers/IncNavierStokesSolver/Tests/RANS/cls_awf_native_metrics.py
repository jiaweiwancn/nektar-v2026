# SPDX-License-Identifier: MIT
"""Shared physical partitions and unambiguous native diagnostic names."""
import math
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import require


def interface_metrics(values):
    require(len(values) == 8, 'unexpected native interface schema')
    return dict(area=values[0], pressure_jump_integral=values[1], normal_traction_jump_RMS=values[2],
        full_vector_traction_jump_RMS=values[3], maximum_absolute_sampled_normal_jump=values[4],
        k_flux_jump_RMS=values[5], normal_velocity_integral=values[6], normal_traction_jump_integral=values[7])


def flat_partition(xml):
    g = ET.parse(xml).getroot().find('GEOMETRY')
    require(g.get('DIM') == '2', 'flat 2D metric only')
    points = [list(map(float, v.text.split())) for v in g.findall('VERTEX/V')]
    require(all(p[2] == 0 for p in points), 'native x/y affine flat chart required')
    return sorted(set(p[0] for p in points)), sorted(set(p[1] for p in points))


def gl4_differences(a, b, xknots):
    """Independent tensor GL4 contraction of piecewise cubic differences.

    Fit actual native GL8 x GL16 samples within EACH reporting panel, verify
    the polynomial representation there, and integrate at GL4 x GL4. No
    essential/matching coefficient is omitted and no field is projected back.
    """
    import numpy as np
    from numpy.polynomial.legendre import leggauss, legvander
    za = a['common_reporting_knots']; require(za == b['common_reporting_knots'], 'reporting partitions differ')
    sa, sb = a['common_samples'], b['common_samples']
    x8, _ = leggauss(8); z16, _ = leggauss(16); x4, w4 = leggauss(4)
    vx, vz, v4 = legvander(x8, 3), legvander(z16, 3), legvander(x4, 3)
    ix, iz = np.linalg.pinv(vx), np.linalg.pinv(vz)
    out = {}; fit_error = {}; area = {}
    for field in ('u', 'v', 'k', 'epsilon', 'pressure'):
        knots = [z for z in za if z >= a['H']] if field == 'epsilon' else za
        nx, nz = len(xknots)-1, len(knots)-1
        delta = (np.array(sa[field])-np.array(sb[field])).reshape(nx, 8, nz, 16)
        blocks = []; weights = []; error = 0.
        for i in range(nx):
            for j in range(nz):
                raw = delta[i, :, j, :]; coeff = ix@raw@iz.T
                error = max(error, float(np.max(abs(vx@coeff@vz.T-raw))))
                blocks.append(v4@coeff@v4.T)
                weights.append(np.outer(w4,w4)*(xknots[i+1]-xknots[i])*(knots[j+1]-knots[j])/4)
        require(error < 2e-11, 'represented difference not cubic on declared panels: '+field)
        values = np.array(blocks); rule = np.array(weights); volume = float(rule.sum())
        if field == 'pressure': values -= float(np.sum(values*rule))/volume
        out[field] = math.sqrt(max(0., float(np.sum(rule*values*values))))
        fit_error[field] = error; area[field] = volume
    return dict(dimensional=out, polynomial_fit_maximum_error=fit_error, physical_areas=area,
        method='independent element/panel GL4 tensor contraction; all represented field values including essential/matching data; one pressure gauge')


def difference_localization(a,b,xknots):
    """Native-panel polynomial reconstruction of K and strong-divergence deltas.

    Derivatives are taken separately in each native physical panel. These are
    GL4 extrema witnesses and L2 norms, not a global polynomial-max certificate.
    """
    import numpy as np
    from numpy.polynomial.legendre import leggauss,legvander,legder,legval
    knots=a['common_reporting_knots'];require(knots==b['common_reporting_knots'],'panels differ')
    x8,_=leggauss(8);z16,_=leggauss(16);g,w=leggauss(4)
    ix=np.linalg.pinv(legvander(x8,3));iz=np.linalg.pinv(legvander(z16,3));v=legvander(g,3)
    dv=np.column_stack([legval(g,legder([0.]*j+[1.])) for j in range(4)])
    nx,nz=len(xknots)-1,len(knots)-1
    delta={key:(np.array(a['common_samples'][key])-np.array(b['common_samples'][key])).reshape(nx,8,nz,16) for key in ('u','v','k')}
    cells=[];divergence=[0.,0.];max_k=(0.,None)
    for i in range(nx):
        for j in range(nz):
            dx=xknots[i+1]-xknots[i];dz=knots[j+1]-knots[j]
            c={key:ix@raw[i,:,j,:]@iz.T for key,raw in delta.items()}
            kval=v@c['k']@v.T;div=2/dx*dv@c['u']@v.T+2/dz*v@c['v']@dv.T
            weights=np.outer(w,w)*dx*dz/4;ks=float(np.sum(weights*kval*kval));ds=float(np.sum(weights*div*div))
            region=0 if knots[j+1]<=a['H'] else 1;divergence[region]+=ds
            p=np.unravel_index(np.argmax(abs(kval)),kval.shape)
            witness=[xknots[i]+(g[p[0]]+1)*dx/2,knots[j]+(g[p[1]]+1)*dz/2]
            if abs(kval[p])>max_k[0]:max_k=(float(abs(kval[p])),witness)
            cells.append(dict(x_bounds=xknots[i:i+2],z_bounds=knots[j:j+2],k_difference_L2_squared=ks,
                strong_divergence_difference_L2_squared=ds))
    return dict(strong_divergence_difference_L2_inner_outer=[math.sqrt(max(0.,x)) for x in divergence],
        largest_k_difference_panels=sorted(cells,key=lambda x:-x['k_difference_L2_squared'])[:8],
        largest_divergence_difference_panels=sorted(cells,key=lambda x:-x['strong_divergence_difference_L2_squared'])[:8],
        sampled_maximum_absolute_k_difference=max_k[0],sampled_k_witness_xz=max_k[1],
        scope='exact represented polynomial derivatives per native/reporting panel; GL4 physical norms; witness is sampled, not an everywhere extremum')
