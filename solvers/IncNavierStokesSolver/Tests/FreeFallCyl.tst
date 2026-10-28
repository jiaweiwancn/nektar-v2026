<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of a free falling circular cylinder</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>FreeFallCyl.xml FreeFallCylc.xml --set-start-time 0 --set-start-chknumber 0</parameters>
    <files>
        <file description="Mesh File">FreeFallCyl.xml</file>
        <file description="Session File">FreeFallCylc.xml</file>
        <file description="initial condition File">FreeFallCyl.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">0.304817</value>
            <value variable="v" tolerance="5e-4">0.303689</value>
            <value variable="p" tolerance="5e-2">0.0166721</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-3">0.151083</value>
            <value variable="v" tolerance="5e-4">0.0797194</value>
            <value variable="p" tolerance="5e-2">0.00747142</value>
        </metric>
    </metrics>
</test>
