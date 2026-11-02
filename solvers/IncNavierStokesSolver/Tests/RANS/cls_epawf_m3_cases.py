"""Reproducible small M3 sessions; generated fields/logs belong on Linux.

Manufactured sources follow the strong conservative equations independently of
the C++ weak operator. The construction has nonzero matching normal velocity,
inner storage/transport, shared full traction and k flux, and outer-only epsilon.
"""
from pathlib import Path
import xml.etree.ElementTree as ET


def set_parameter(root, name, value):
    parameters = root.find("CONDITIONS/PARAMETERS")
    for p in parameters:
        if p.text.split("=", 1)[0].strip() == name:
            p.text = f"{name} = {value}"
            return
    ET.SubElement(parameters, "P").text = f"{name} = {value}"


def set_info(root, name, value):
    info = root.find("CONDITIONS/SOLVERINFO")
    for p in info:
        if p.attrib["PROPERTY"] == name:
            p.set("VALUE", str(value))
            return
    ET.SubElement(info, "I", {"PROPERTY": name, "VALUE": str(value)})


def mesh(root, knots, degree, points=None):
    points = points or max(8, degree+3)
    for e in root.find("EXPANSIONS"):
        p = degree-2 if e.attrib["FIELDS"] == "p" else degree
        e.set("NUMMODES", f"{p+1},{p+1}")
        e.set("NUMPOINTS", f"{points},{points}")
    geometry = root.find("GEOMETRY")
    geometry.clear()
    geometry.attrib.update({"DIM": "2", "SPACE": "2"})
    vertices, edges, elements = (ET.SubElement(geometry, name) for name in ("VERTEX", "EDGE", "ELEMENT"))
    for j, r in enumerate(knots):
        for i in (0, 1):
            ET.SubElement(vertices, "V", {"ID": str(2*j+i)}).text = f"{i} {r:.17g} 0"
    edge_ids = {}
    for j in range(len(knots)-1):
        v = [2*j, 2*j+1, 2*j+3, 2*j+2]
        ids = []
        for a, b in zip(v, v[1:]+v[:1]):
            key = tuple(sorted((a, b)))
            if key not in edge_ids:
                edge_ids[key] = len(edge_ids)
                ET.SubElement(edges, "E", {"ID": str(edge_ids[key])}).text = f"{a} {b}"
            ids.append(edge_ids[key])
        ET.SubElement(elements, "Q", {"ID": str(j)}).text = " ".join(map(str, ids))
    composite = ET.SubElement(geometry, "COMPOSITE")
    ET.SubElement(composite, "C", {"ID": "0"}).text = f"Q[0-{len(knots)-2}]"
    ET.SubElement(geometry, "DOMAIN").text = "C[0]"


def write(root, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(destination, encoding="utf-8", xml_declaration=True)
    return destination


def parallel_session(template, destination, degree, inner_knots, normal_order=16):
    root = ET.parse(template).getroot()
    mesh(root, list(inner_knots)+[4.], degree)
    set_parameter(root, "RANSWallNormalOrder", normal_order)
    return write(root, destination)


def evolving_session(template, destination, degree=6, time=.3, mode="Solve", condensed=False,
                     alpha=5., perturb=True):
    import sympy as sp
    from sympy.printing.str import StrPrinter

    class NektarPrinter(StrPrinter):
        def _print_Pow(self, expression):
            base, exponent = expression.as_base_exp()
            return f"({self._print(base)})^({self._print(exponent)})"
        def _print_Abs(self, expression):
            return "fabs("+self._print(expression.args[0])+")"

    printer = NektarPrinter()
    root = ET.parse(template).getroot()
    # A thin first outer element actually samples the retained active Yap branch.
    mesh(root, [0., 2., 2.02, 4.], degree)
    set_parameter(root, "Time", repr(time))
    set_parameter(root, "RANSWallStageAlpha", repr(alpha))
    set_parameter(root, "RANSWallTangentialOrder", max(8, degree+3))
    set_parameter(root, "M3Top", 4)
    set_info(root, "RANSWallCondense", "True" if condensed else "False")
    set_info(root, "RANSWallStageMode", mode)
    conditions = root.find("CONDITIONS")
    for old in list(conditions.findall("FUNCTION")):
        conditions.remove(old)
    x, r, t = sp.symbols("x y t", real=True)
    H, Z, rK, rough = sp.symbols("RANSWallHeight M3Top RANSWallSampleHeight RANSWallRoughness", positive=True)
    K0, E0, slope, nu, qstar, stage_alpha = sp.symbols("K0 E0 Slope Kinvis Qstar RANSWallStageAlpha", positive=True)
    zv = sp.Symbol("Zv", real=True)
    s = r/H
    gate = 6*s**2-8*s**3+3*s**4
    primitive = r*r/H-r**3/(3*H**2)
    U = sp.Rational(1, 5)+sp.Rational(3, 100)*sp.sin(t)
    W = sp.Rational(1, 100)+sp.Rational(1, 200)*sp.cos(t)
    u = -W*x*sp.diff(gate, r)+U*sp.diff(primitive, r)
    w = W*gate
    # The side traces are exactly K0 at every p. This keeps the independently
    # projected outer epsilon corner compatible with its sampled-K trace even
    # in the deliberately unresolved p=4 audit.
    ki = K0+4*x*(1-x)*(sp.Rational(1, 25)+sp.sin(t)/100)*r**2*(r-rK)**2*(r-H)**2/H**6
    ko = K0+(sp.Rational(1, 50)+sp.cos(t)/200)*(r-H)**2/(Z-H)**2
    eo = E0+(sp.Rational(3, 1000)+sp.sin(t)/1000)*(r-H)**2/(Z-H)**2
    px, pzi, pzo = sp.Rational(1, 50)*(1+sp.sin(t)/10), sp.Rational(3, 100)*sp.cos(t), sp.sin(t)/25
    mean = (pzi*H**3+pzo*(Z-H)**3)/(3*Z)
    pi, po = px*(x-sp.Rational(1, 2))+pzi*(r-H)**2-mean, px*(x-sp.Rational(1, 2))+pzo*(r-H)**2-mean
    nt, D = slope*(r-zv), nu+slope*(r-zv)
    production = nt*(2*sp.diff(u, x)**2+2*sp.diff(w, r)**2+sp.diff(u, r)**2)
    p0 = nt*(K0/qstar)**2/D**2
    # F_A = px - b_outer,x = -dU_H/dt. The roughness force is
    # b_inner + (1-a_h)F_A; no mass or residual is multiplied by a_h.
    below = (1-sp.sign(r-rough))/2
    bi_u = sp.diff(u, t)+u*sp.diff(u, x)+w*sp.diff(u, r)+sp.diff(pi, x)-sp.diff(D*sp.diff(u, r), r)+below*sp.diff(U, t)
    bi_w = sp.diff(w, t)+w*sp.diff(w, r)+sp.diff(pi, r)-D*sp.diff(w, r, 2)-2*slope*sp.diff(w, r)
    bo_u, bo_w = sp.diff(U, t)+sp.diff(po, x), sp.diff(W, t)+sp.diff(po, r)
    source_ki = (sp.diff(ki, t)+u*sp.diff(ki, x)+w*sp.diff(ki, r)
                 -D*sp.diff(ki, x, 2)-sp.diff(D*sp.diff(ki, r), r)-production+p0)
    cmu0 = sp.Rational(3, 10)*(1-sp.exp(-sp.Rational(9, 25)))
    nto = cmu0*ko**2/eo
    source_ko = sp.diff(ko, t)+W*sp.diff(ko, r)-sp.diff(nto*sp.diff(ko, r), r)+eo
    ratio = ko**sp.Rational(3, 2)/(sp.Rational(5, 2)*eo*r)
    positive = ((ratio-1)+sp.Abs(ratio-1))/2
    yap = sp.Rational(83, 100)*eo**2/ko*positive*ratio**2
    source_eo = sp.diff(eo, t)+W*sp.diff(eo, r)-sp.diff(nto*sp.diff(eo, r)/sp.Rational(13, 10), r)+sp.Rational(48, 25)*eo**2/ko-yap
    exact_i, exact_o = {"u": u, "v": w, "k": ki, "p": pi}, {"u": U, "v": W, "k": ko, "eps_tilde": eo, "p": po}
    def function(name, values):
        node = ET.SubElement(conditions, "FUNCTION", {"NAME": name})
        for variable, expression in values.items():
            ET.SubElement(node, "E", {"VAR": variable, "VALUE": printer.doprint(expression)})
    for side, exact in (("Inner", exact_i), ("Outer", exact_o)):
        initial = dict(exact)
        if perturb and mode == "Solve":
            psi = sp.Rational(1, 100000)*x**2*(1-x)**2*r**2*(Z-r)**2
            initial["u"] += sp.diff(psi, r)
            initial["v"] -= sp.diff(psi, x)
        function("InitialConditions"+side, initial)
        function("ExactSolution"+side, exact)
        function("BoundaryValues"+side, {v: e for v, e in exact.items() if v != "p"})
        function("StageHistory"+side, {v: e-sp.diff(e, t)/stage_alpha if v != "p" else e for v, e in exact.items()})
    function("BodyForceInner", {"u": bi_u, "v": bi_w})
    function("BodyForceOuter", {"u": bo_u, "v": bo_w})
    function("ScalarSourceInner", {"k": source_ki})
    function("ScalarSourceOuter", {"k": source_ko, "eps_tilde": source_eo})
    return write(root, destination)
