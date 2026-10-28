<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Full 3D HIT inlet, zero RMS, periodic transverse sides and HOutflow</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>m3-mean-only.xml</parameters>
    <files>
        <file description="Session File">HITInflow/m3-mean-only.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">0</value>
            <value variable="v" tolerance="1e-10">0</value>
            <value variable="w" tolerance="1e-10">0</value>
            <value variable="p" tolerance="1e-10">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">0</value>
            <value variable="v" tolerance="1e-10">0</value>
            <value variable="w" tolerance="1e-10">0</value>
            <value variable="p" tolerance="1e-10">0</value>
        </metric>
    </metrics>
</test>
