<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>HIT inlet restart through Checkpoint filter, compared with uninterrupted IMEX flow</description>
    <segment type="sequential">
        <executable>IncNavierStokesSolver</executable>
        <parameters>m4-h1d-von-karman.xml</parameters>
        <processes>1</processes>
    </segment>
    <segment type="sequential">
        <executable>IncNavierStokesSolver</executable>
        <parameters>m4-h1d-checkpoint-filter.xml</parameters>
        <processes>1</processes>
    </segment>
    <segment type="sequential">
        <executable>IncNavierStokesSolver</executable>
        <parameters>m4-h1d-filter-restart.xml</parameters>
        <processes>1</processes>
    </segment>
    <segment type="sequential">
        <executable>HITFieldCompare</executable>
        <parameters>m4-h1d-von-karman.xml m4-h1d-von-karman.fld m4-h1d-filter-restart.fld</parameters>
        <processes>1</processes>
    </segment>
    <files>
        <file description="HIT regression input">HITInflow/m4-h1d-von-karman.xml</file>
        <file description="HIT regression input">HITInflow/m4-h1d-checkpoint-filter.xml</file>
        <file description="HIT regression input">HITInflow/m4-h1d-filter-restart.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>"([uvwp])":\{"l2":[^,]+,"linf":([^,]+),.*</regex>
            <matches>
                <match>
                    <field id="0">u</field>
                    <field id="1" tolerance="0.0002">0</field>
                </match>
                <match>
                    <field id="0">v</field>
                    <field id="1" tolerance="0.0002">0</field>
                </match>
                <match>
                    <field id="0">w</field>
                    <field id="1" tolerance="0.0002">0</field>
                </match>
                <match>
                    <field id="0">p</field>
                    <field id="1" tolerance="0.002">0</field>
                </match>
            </matches>
        </metric>
        <metric type="Regex" id="2">
            <regex>HIT_FIELD_IDENTITY identical ([01]) same_time ([01])</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                    <field id="1">1</field>
                </match>
            </matches>
        </metric>
        <metric type="Regex" id="3">
            <regex>HIT_FIELD_NONZERO velocity_components ([0-3])</regex>
            <matches>
                <match>
                    <field id="0">3</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
