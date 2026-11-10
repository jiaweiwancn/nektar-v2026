#!/usr/bin/env python3
import copy,json,tempfile,unittest
from pathlib import Path
from rans_request_manifest import RequestManifest
class ManifestTests(unittest.TestCase):
 def test_reject_before_mutating(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)/'job';request=dict(quadrature=[32,8],perturbation=0.,horizon=.016,runtime='binary-A',boundary='map-A',evaluation=dict(sides=[-1,1]))
   m=RequestManifest(root,request);(root/'results.json').write_text('unchanged');before={p.name:p.read_bytes() for p in root.iterdir()}
   with self.assertRaises(FileExistsError):RequestManifest(root,request)
   for key,value in [('quadrature',[64,16]),('perturbation',.01),('horizon',.1),('runtime','binary-B'),('boundary','map-B'),('evaluation',dict(sides=[1,-1]))]:
    other=copy.deepcopy(request);other[key]=value
    with self.assertRaisesRegex(ValueError,'request mismatch'):RequestManifest(root,other,True)
    self.assertEqual(before,{p.name:p.read_bytes() for p in root.iterdir()})
 def test_identical_resume_missing_only(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)/'job';request=dict(physics='E',runtime='immutable');m=RequestManifest(root,request)
   first=root/'first';first.mkdir();(first/'state.bin').write_bytes(b'preserved state');m.seal(first,['state.bin']);snapshot=(first/'state.bin').stat().st_mtime_ns
   resumed=RequestManifest(root,request,True);self.assertTrue(resumed.completed(first));self.assertFalse(resumed.completed(root/'second'))
   second=root/'second';second.mkdir();(second/'state.bin').write_bytes(b'new work');resumed.seal(second,['state.bin'])
   self.assertEqual(snapshot,(first/'state.bin').stat().st_mtime_ns)
   (first/'state.bin').write_bytes(b'corrupt')
   with self.assertRaisesRegex(ValueError,'evidence changed'):resumed.completed(first)
if __name__=='__main__':unittest.main()
