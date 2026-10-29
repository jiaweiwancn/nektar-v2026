<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Manufactured HIT inlet: evolved bulk divergence, net flux and Gauss consistency</description>
    <segment type="sequential">
        <executable>IncNavierStokesSolver</executable>
        <parameters>m3-manufactured-dt0.01000.xml</parameters>
        <processes>1</processes>
    </segment>
    <segment type="sequential">
        <executable>HITFlowDiagnostics</executable>
        <parameters>m3-manufactured-dt0.01000.xml m3-manufactured-dt0.01000.fld</parameters>
        <processes>1</processes>
    </segment>
    <files>
        <file description="Manufactured HIT input">HITInflow/m3-manufactured-dt0.01000.xml</file>
        <file description="Manufactured HIT input">HITInflow/m3-narrow-band.dat</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>HIT_FLOW div_l2 ([^ ]+) div_relative ([^ ]+) flux_relative ([^ ]+) gauss_relative ([^ ]+)</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-05">0</field>
                    <field id="1" tolerance="1e-05">0</field>
                    <field id="2" tolerance="1e-08">0</field>
                    <field id="3" tolerance="1e-11">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
