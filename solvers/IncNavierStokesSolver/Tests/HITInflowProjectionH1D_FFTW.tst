<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>HIT constrained projection, nonzero packed Fourier components and repeatability</description>
    <executable>HITBoundaryProbe</executable>
    <parameters>m4-h1d-fftw.xml</parameters>
    <files>
        <file description="HIT regression input">HITInflow/m4-h1d-fftw.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>HIT_PROBE raw_error ([^ ]+) projection_error ([^ ]+) nonzero_z ([01]) repeat ([01])</regex>
            <matches>
                <match>
                    <field id="0" tolerance="5e-13">0</field>
                    <field id="1" tolerance="2e-06">0</field>
                    <field id="2">1</field>
                    <field id="3">1</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
