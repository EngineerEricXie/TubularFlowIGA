"""Exercise rank-layout reassembly independently of solver output."""
import math
from pathlib import Path
import struct
import tempfile
import unittest

from hpc_compare_checkpoint_fields import owned_field
from hpc_solver_prefixes import norms


class OwnedFields(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def shard(self, name, total, begin, end, values):
        path = self.root/name
        path.write_bytes(b'header\n'*6+struct.pack('<QQQ', total, begin, end)+struct.pack('<'+'d'*len(values), *values))
        return path

    def test_reordered_partition_and_empty_rank(self):
        paths = [self.shard('last', 4, 2, 4, [3.,4.]), self.shard('first', 4, 0, 2, [1.,2.]),
                 self.shard('empty', 4, 4, 4, [])]
        self.assertEqual(owned_field(paths), [1.,2.,3.,4.])

    def test_reject_invalid_partition(self):
        first = self.shard('first', 4, 0, 2, [1.,2.])
        for total,begin,end,values in [(4,1,4,[2.,3.,4.]), (4,3,4,[4.]), (5,2,4,[3.,4.]),
                                       (4,2,5,[3.,4.,5.]), (4,2,4,[3.])]:
            with self.subTest(total=total, begin=begin, end=end, values=values):
                with self.assertRaises(ValueError):
                    owned_field([first, self.shard('other', total, begin, end, values)])
        with self.assertRaises(ValueError):
            owned_field([first])
        with self.assertRaises(ValueError):
            owned_field([])

    def test_nonfinite_and_wrong_values_fail_gate(self):
        for values in ([math.nan], [math.inf], [2.]):
            with self.assertRaises(AssertionError):
                norms([1.], values)


if __name__ == '__main__':
    unittest.main()
