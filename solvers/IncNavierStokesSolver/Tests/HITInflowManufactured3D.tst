<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Full 3D seeded time-dependent HIT inlet with manufactured body force</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>m3-manufactured-dt0.01000.xml</parameters>
    <files>
        <file description="Session File">HITInflow/m3-manufactured-dt0.01000.xml</file>
        <file description="Shell Spectrum">HITInflow/m3-narrow-band.dat</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0</value>
            <value variable="v" tolerance="3e-5">0</value>
            <value variable="w" tolerance="2e-5">0</value>
            <value variable="p" tolerance="1e-5">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">0</value>
            <value variable="v" tolerance="5e-5">0</value>
            <value variable="w" tolerance="1e-5">0</value>
            <value variable="p" tolerance="1e-5">0</value>
        </metric>
    </metrics>
</test>
