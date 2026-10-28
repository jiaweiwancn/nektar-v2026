<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with rotation, simple domain without body</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>MovingRefFrame_Rot_SimpleDomain.xml</parameters>
    <files>
        <file description="Session File">MovingRefFrame_Rot_SimpleDomain.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">0</value>
            <value variable="v" tolerance="5e-4">0</value>
            <value variable="p" tolerance="5e-2">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-3">0</value>
            <value variable="v" tolerance="5e-4">0</value>
            <value variable="p" tolerance="2e-2">0</value>
        </metric>
    </metrics>
</test>
