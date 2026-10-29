<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Seeded homogeneous HIT inlet with tabulated shell spectrum</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>m4-h1d-tabulated.xml</parameters>
    <files>
        <file description="HIT regression input">HITInflow/m4-h1d-tabulated.xml</file>
        <file description="HIT regression input">HITInflow/m4-shell-spectrum.dat</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="2e-5">6.28855</value>
            <value variable="v" tolerance="2e-6">0.139403</value>
            <value variable="w" tolerance="2e-6">0.155972</value>
            <value variable="p" tolerance="2e-6">0.00803102</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-5">1.10136</value>
            <value variable="v" tolerance="2e-6">0.0986928</value>
            <value variable="w" tolerance="2e-6">0.108804</value>
            <value variable="p" tolerance="2e-6">0.00675348</value>
        </metric>
    </metrics>
</test>
