# SPDX-License-Identifier: MIT
import tempfile,unittest
from pathlib import Path
from wj_h15_observer import atomic,identity
from wj_h18_readback import selection

class VerifiedReadback(unittest.TestCase):
    def fixture(self,r,name,step,startup=False):
        m=r/f'science/target/{name}/monitor';m.mkdir(parents=True)
        files={}
        for k in ('checkpoint','callback','report'):
            p=m/k;p.write_bytes(k.encode());files[k]=identity(p)
        if startup:
            p=m/'completion.json';atomic(p,files)
            q=dict(schema='WJStartupVerificationCommitV1',binding={'absolute_step':str(step)},completion=identity(p))
            commit=m/'startup.verification-commit.json'
        else:
            q=dict(schema='WJQualificationCommitV1',binding={'absolute_step':str(step)},**files);commit=m/f'step-{step}.qualification-commit.json'
        atomic(commit,q);atomic(m/'latest-fully-verified-checkpoint.json',q)
        return m

    def test_restarted_input_is_last_verified_not_earlier_branch_future(self):
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);self.fixture(r,'A',4);self.fixture(r,'B',2,True)
            folder,request=selection(r)
            self.assertEqual(folder.name,'B');self.assertTrue(request['qualified_cursor']['path'].endswith('startup.verification-commit.json'))
            self.assertIn('/B/',request['native_report'])

    def test_prepared_pointer_cannot_authorize_and_corrupt_pin_is_not_fallback(self):
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);m=self.fixture(r,'A',4)
            atomic(m/'latest-fully-verified-checkpoint.json',{'schema':'PREPARED_PASS'},False)
            with self.assertRaises(AssertionError):selection(r)
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);m=self.fixture(r,'B',2,True);(m/'callback').write_bytes(b'changed')
            with self.assertRaises(AssertionError):selection(r)

    def test_completed_commit_and_original_cp0_roles(self):
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);self.fixture(r,'A',4)
            self.assertEqual(selection(r)[0].name,'A')
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);self.fixture(r,'A',0,True);self.assertIsNone(selection(r))

if __name__=='__main__':unittest.main()
