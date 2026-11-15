# SPDX-License-Identifier: MIT
import tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from publish_wj_h20_resumed import saved_arrays

class Publication(unittest.TestCase):
    def test_saved_array_reads_are_once_per_key_and_preserve_IEEE_bytes(self):
        arrays=dict(U=np.array([[0.,-0.,1.],[2.,3.,4.]]),geometry=np.array([251,265],dtype=np.int64))
        calls=[]
        class Archive:
            files=list(arrays)
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def __getitem__(self,key):calls.append(key);return arrays[key]
        with patch('publish_wj_h20_resumed.np.load',return_value=Archive()):loaded=saved_arrays('immutable.npz')
        for _ in range(10):self.assertEqual(loaded['U'].tobytes(),arrays['U'].tobytes())
        self.assertEqual(calls,list(arrays))
        with tempfile.TemporaryDirectory()as t:
            p=Path(t)/'actual.npz';np.savez_compressed(p,**arrays);loaded=saved_arrays(p)
            for key in arrays:self.assertEqual(loaded[key].tobytes(),arrays[key].tobytes())

if __name__=='__main__':unittest.main()
