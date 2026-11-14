# SPDX-License-Identifier: MIT
import copy
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch
from wj_h18_observer import due,rotate_recoveries,rotate_owned,parity
from study_wj_h18 import protected,xml_case
from run_wj_h18 import small_selection,CAPS,ROOT
from wj_h17_observer import equal

class CheckpointAudits(unittest.TestCase):
    def test_absolute_cadence_restart_final_and_stride_one(self):
        i=dict(scientific_stride=100,maximum_absolute_step=202010,terminal_start=180000,terminal_stride=50)
        self.assertFalse(due(99,i));self.assertTrue(due(100,i));self.assertTrue(due(180050,i));self.assertTrue(due(202010,i))
        self.assertEqual([s for s in range(3,5)if due(s,dict(scientific_stride=4,maximum_absolute_step=4))],[4])
        self.assertTrue(all(due(s,dict(scientific_stride=1,maximum_absolute_step=100))for s in range(1,101)))

    def test_recovery_rotation_preserves_science_pending_failure_and_verified(self):
        rows=[]
        protected_ids={'science','pending','failure','verified'}
        for i in range(10000):
            rows.append(dict(id=str(i),step=i,role='recovery'))
            removed=rotate_recoveries(rows,protected_ids,3);ids={r['id']for r in removed};rows=[r for r in rows if r['id']not in ids]
            self.assertLessEqual(len(rows),3)
        for i,key in enumerate(protected_ids):rows.append(dict(id=key,step=i,role='recovery'))
        self.assertTrue(all(v['id']not in protected_ids for v in rotate_recoveries(rows,protected_ids,3)))
        rows.append(dict(id='archived-scientific',step=-1,role='scientific'))
        self.assertNotIn('archived-scientific',{v['id']for v in rotate_recoveries(rows,set(),3)})

    def test_owned_rotation_is_durable_and_never_deletes_pending(self):
        from wj_h15_observer import identity
        with tempfile.TemporaryDirectory()as d:
            root=Path(d);records=[]
            for n in range(12):
                f=root/str(n);f.write_bytes(bytes([n]));records.append(dict(id=str(n),step=n,role='recovery',verification='NOT_REQUIRED_SUPERSEDED',payloads=[identity(f)]))
            removed=rotate_owned(root,records,{'0'},3)
            self.assertEqual(removed,[str(n)for n in range(1,9)])
            self.assertTrue((root/'0').exists());self.assertTrue(all((root/str(n)).exists()for n in (9,10,11)))
            records=[records[0]]+records[9:];records[0]['verification']='PENDING'
            with self.assertRaisesRegex(AssertionError,'pending'):rotate_owned(root,records,set(),2)

    def test_five_successful_and_one_deliberate_failure_reservations(self):
        import sys
        root=Path('/home/jay/build')/ROOT
        counts={'legacy':2,'A':2,'B':1,'fault':1}
        for name,count in counts.items():
            cmd=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h18.py')),'--root',str(root),'--mode','small','--case',name]
            self.assertEqual(small_selection(root,cmd,count),name)
            with self.assertRaises(AssertionError):small_selection(root,cmd,count+1)
        self.assertEqual(sum(counts[n]for n in ('legacy','A','B')),5)
        self.assertEqual(CAPS['build'],1200);self.assertEqual(CAPS['qualification'],3000)

    def test_physical_xml_time_functions_and_dt_are_not_execution_controls(self):
        x=ET.fromstring('<NEKTAR><CONDITIONS><PARAMETERS><P>Time = 0</P><P>TimeStep = .0025</P><P>NumSteps = 4</P></PARAMETERS><SOLVERINFO><I PROPERTY="RANSWallMonitoringPolicy" VALUE="Off"/></SOLVERINFO><FUNCTION NAME="InitialConditions"><E VALUE="log(y)"/></FUNCTION></CONDITIONS></NEKTAR>')
        y=copy.deepcopy(x);y.find('./CONDITIONS/PARAMETERS/P[3]').text='NumSteps = 2';y.find('./CONDITIONS/SOLVERINFO/I').set('VALUE','CheckpointAuditV1')
        self.assertEqual(protected(x),protected(y))
        for node,text in (('./CONDITIONS/PARAMETERS/P[1]','Time = .005'),('./CONDITIONS/PARAMETERS/P[2]','TimeStep = .005')):
            z=copy.deepcopy(y);z.find(node).text=text;self.assertNotEqual(protected(x),protected(z))
        z=copy.deepcopy(y);z.find('./CONDITIONS/FUNCTION/E').set('VALUE','sqrt(y)');self.assertNotEqual(protected(x),protected(z))

    def test_signed_zero_is_part_of_physical_parity(self):
        with self.assertRaises(AssertionError):equal([0.],[-0.],'serialized coefficients')

    def test_pending_is_distinct_from_qualified_and_prepared_cannot_promote(self):
        from study_wj_h15 import qualified_step
        from wj_h15_observer import atomic,identity
        with tempfile.TemporaryDirectory()as d:
            f=Path(d);m=f/'monitor';m.mkdir();p=m/'payload';p.write_bytes(b'original')
            atomic(m/'latest-durable-checkpoint.json',dict(status='NATIVE_COMPLETED_PENDING_VERIFICATION',absolute_step=2))
            atomic(m/'step-2.ready.qualified.json',dict(status='PREPARED_PASS'))
            self.assertEqual(qualified_step(f),0)

if __name__=='__main__':unittest.main()
