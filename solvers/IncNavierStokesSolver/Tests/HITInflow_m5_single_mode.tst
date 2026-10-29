<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>One seeded transverse Fourier pair with exact advected velocity and body force; manufactured accuracy test</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>m5-single-mode.xml</parameters>
    <files>
        <file description="Session File">HITInflow/m5-single-mode.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">0</value>
            <value variable="v" tolerance="2e-6">0</value>
            <value variable="w" tolerance="2e-6">0</value>
            <value variable="p" tolerance="1e-9">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">0</value>
            <value variable="v" tolerance="3e-6">0</value>
            <value variable="w" tolerance="3e-6">0</value>
            <value variable="p" tolerance="1e-9">0</value>
        </metric>
        <metric type="regex" id="3">
            <regex>^HITInflow independent modes: ([0-9]+),.*$</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
