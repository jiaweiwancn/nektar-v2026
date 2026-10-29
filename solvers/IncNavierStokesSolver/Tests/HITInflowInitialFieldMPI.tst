<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Native HIT initial-field projection preserves a constant mean through mesh/Fourier MPI assembly</description>
    <segment type="sequential">
        <executable>HITInitialField</executable>
        <parameters>--npz 2 m4-h1d-mean-only.xml initial.fld</parameters>
        <processes>4</processes>
    </segment>
    <segment type="sequential">
        <executable>IncNavierStokesSolver</executable>
        <parameters>--npz 2 m4-h1d-mean-only.xml</parameters>
        <processes>4</processes>
    </segment>
    <segment type="sequential">
        <executable>HITFieldCompare</executable>
        <parameters>m4-h1d-mean-only.xml initial.fld m4-h1d-mean-only.fld</parameters>
        <processes>1</processes>
    </segment>
    <files>
        <file description="Mean initial-field assembly input">HITInflow/m4-h1d-mean-only.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>"([uvwp])":\{"l2":[^,]+,"linf":([^,]+),.*</regex>
            <matches>
                <match><field id="0">u</field><field id="1" tolerance="5e-12">0</field></match>
                <match><field id="0">v</field><field id="1" tolerance="5e-12">0</field></match>
                <match><field id="0">w</field><field id="1" tolerance="5e-12">0</field></match>
                <match><field id="0">p</field><field id="1" tolerance="5e-12">0</field></match>
            </matches>
        </metric>
    </metrics>
</test>
