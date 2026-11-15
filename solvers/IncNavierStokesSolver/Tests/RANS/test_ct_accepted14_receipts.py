# SPDX-License-Identifier: MIT
"""Discriminate forged capacity/accuracy labels and missing actual solve bindings."""
import copy,json,tempfile,unittest
from pathlib import Path
from ct_accepted14_pair import numerical_receipts

class ReceiptTests(unittest.TestCase):
    def setUp(self):
        self.valid=[dict(event='factor',id=1,rows=48,matrix_nnz=140,factor_nnz=300,mass=False,ordering='VectorCOLAMDV1'),dict(event='solve',id=1,rhs_maximum=1e-18,true_residual=1e-30,bound=2e-30+1e-25,refinements=3)]
    def evaluate(self,rows):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'receipts.jsonl';p.write_text(''.join(json.dumps(x)+'\n' for x in rows));return numerical_receipts(p)
    def test_tiny_actual_rhs(self):
        self.assertEqual(self.evaluate(self.valid)['status'],'pass')
    def test_rejects_invalid_actual_primitives(self):
        for label,index,key,value in [('fill',0,'factor_nnz',40000001),('nnz',0,'matrix_nnz',20000001),('foreign_factor',1,'id',2),('residual',1,'true_residual',1e-20),('too_many_refinements',1,'refinements',4),('nonfinite',1,'rhs_maximum',float('nan'))]:
            with self.subTest(label=label):
                rows=copy.deepcopy(self.valid);rows[index][key]=value
                with self.assertRaises(ValueError):self.evaluate(rows)
    def test_rejects_missing_or_duplicate_evidence(self):
        for rows in [[],self.valid[:1],self.valid[1:],[self.valid[0],*self.valid]]:
            with self.subTest(rows=rows):
                with self.assertRaises(ValueError):self.evaluate(rows)

if __name__=='__main__':unittest.main()
