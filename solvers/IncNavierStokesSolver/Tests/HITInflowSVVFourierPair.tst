<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Homogeneous Fourier SVV preserves cosine/sine partner damping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>m6-svv-fourier-pair.xml</parameters>
    <files>
        <file description="Session File">HITInflow/m6-svv-fourier-pair.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0</value>
            <value variable="v" tolerance="1e-7">0</value>
            <value variable="w" tolerance="1e-7">0</value>
            <value variable="p" tolerance="1e-7">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0</value>
            <value variable="v" tolerance="1e-7">0</value>
            <value variable="w" tolerance="1e-7">0</value>
            <value variable="p" tolerance="1e-7">0</value>
        </metric>
    </metrics>
</test>
