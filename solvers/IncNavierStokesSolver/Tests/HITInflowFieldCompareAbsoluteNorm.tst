<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>HIT field comparison reports absolute error for a uniformly negative difference</description>
    <segment type="sequential">
        <executable>HITInitialField</executable>
        <parameters>--parameter Ubar=1 m6-h1d-mean-initial.xml one.fld</parameters>
        <processes>1</processes>
    </segment>
    <segment type="sequential">
        <executable>HITInitialField</executable>
        <parameters>--parameter Ubar=2 m6-h1d-mean-initial.xml two.fld</parameters>
        <processes>1</processes>
    </segment>
    <segment type="sequential">
        <executable>HITFieldCompare</executable>
        <parameters>m6-h1d-mean-initial.xml one.fld two.fld</parameters>
        <processes>1</processes>
    </segment>
    <files>
        <file description="Mean-only initial field input">HITInflow/m6-h1d-mean-initial.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <!-- Regex numeric tolerances compare magnitudes, so reject a
                 leading minus in the output as well as checking its value. -->
            <regex>"([uvwp])":\{"l2":[^,]+,"linf":([0-9][^,]*),.*</regex>
            <matches>
                <match>
                    <field id="0">u</field>
                    <field id="1" tolerance="5e-12">1</field>
                </match>
                <match>
                    <field id="0">v</field>
                    <field id="1" tolerance="5e-12">0</field>
                </match>
                <match>
                    <field id="0">w</field>
                    <field id="1" tolerance="5e-12">0</field>
                </match>
                <match>
                    <field id="0">p</field>
                    <field id="1" tolerance="5e-12">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
