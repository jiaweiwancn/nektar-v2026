<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Seeded homogeneous HIT inlet with von Karman shell spectrum</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>m4-h1d-von-karman.xml</parameters>
    <files>
        <file description="HIT regression input">HITInflow/m4-h1d-von-karman.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="2e-5">6.28892</value>
            <value variable="v" tolerance="2e-6">0.100824</value>
            <value variable="w" tolerance="2e-6">0.195365</value>
            <value variable="p" tolerance="2e-6">0.00817315</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-5">1.13202</value>
            <value variable="v" tolerance="2e-6">0.059959</value>
            <value variable="w" tolerance="2e-6">0.125703</value>
            <value variable="p" tolerance="2e-6">0.00587705</value>
        </metric>
    </metrics>
</test>
