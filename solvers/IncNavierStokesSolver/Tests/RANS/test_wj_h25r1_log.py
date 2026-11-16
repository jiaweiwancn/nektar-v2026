# SPDX-License-Identifier: MIT
import hashlib,io,tempfile,unittest
from pathlib import Path
from bounded_wj_h25r1_log import BoundedLog


class BoundedStdout(unittest.TestCase):
    def test_rotation_preserves_recent_bytes_prefix_and_full_stream_hash(self):
        data=bytes(range(251))*277
        with tempfile.TemporaryDirectory()as tmp:
            p=Path(tmp)/'native.log';log=BoundedLog(p,chunk=4096,keep=2,prefix=128)
            log.attach(io.BytesIO(data));r=log.finish()
            self.assertEqual(r['full_stream_sha256'],hashlib.sha256(data).hexdigest())
            self.assertEqual(r['total_emitted_bytes'],len(data));self.assertGreater(r['rotations'],2)
            self.assertEqual((p.with_name(p.name+'.first')).read_bytes(),data[:128])
            recent=b''.join(p.with_name(p.name+'.'+str(i)).read_bytes()for i in (2,1))+p.read_bytes()
            self.assertEqual(recent,data[-len(recent):])
            self.assertLessEqual(sum(v['bytes']for v in r['retained']),3*4096+128)

    def test_existing_owned_log_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory()as tmp:
            p=Path(tmp)/'native.log';p.write_bytes(b'preserve')
            with self.assertRaises(AssertionError):BoundedLog(p)
            self.assertEqual(p.read_bytes(),b'preserve')


if __name__=='__main__':unittest.main()
