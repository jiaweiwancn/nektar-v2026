The CSV contains sixty hollow rough-case symbols from Cao–Tamura Fig.3(a–d),
printed/PDF p.7. White-filled vector circle centers were extracted from the
original PDF content stream and visually checked with red cross overlays.
The black filled smooth-case symbols, Raupach crosses and legend symbols are
excluded. `digitization.json` retains labeled tick calibration, source hash,
marker coordinates, normalizations and a conservative reading resolution.

Height origins differ: panel(a) uses physical y in mm; (b) uses y/delta;
(c,d) use (y-d)/delta. Convert with delta=250mm and the paper's assumed d=3.5mm.
The bottom border of panel(b) is not zero. Do not apply its border as a zero
calibration or forget the displacement in(c,d).

These are digitized values, separate from rounded reported near-wall anchors.
The reading resolution is not an experimental confidence interval. The paper's
probe-voltage calibration percentage is not assigned to every velocity sample.
The Fig.3 plane is the virtual crest with the hill absent. No hill measurement
is fitted here, and no spanwise variance or epsilon is asserted as measured.

Reproduce with `digitize_cao_tamura_fig3.py` and the hash-matched original PDF.
The overlay is a checking artifact; the PDF is not copied into this repository.
