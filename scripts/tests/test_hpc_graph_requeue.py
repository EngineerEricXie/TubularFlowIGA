import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[2]


def run_wrapper(root):
    original = (ROOT / "solvers/cpu/slurm/multinode_graph.sbatch").read_text()
    results = {}
    for variant, script in [('actual', original)]:
        directory = root/variant
        directory.mkdir()
        repo, binaries = directory/'repo', directory/'bin'
        (repo/'scripts').mkdir(parents=True)
        (repo/'solvers/coupling').mkdir(parents=True)
        binaries.mkdir()
        (directory/'case').mkdir()
        (directory/'spooled-script').write_text(script)
        dummy = repo/'solvers/coupling/iga_multidomain_flow'
        dummy.write_text('#!/bin/bash\nexit 0\n'); dummy.chmod(0o755)
        recorder = '''
    import argparse,json
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--output');p.add_argument('--status');p.add_argument('--returncode');p.add_argument('--checkpoint-requested',action='store_true')
    a,_=p.parse_known_args();Path(a.output).write_text(json.dumps(vars(a)))
    '''
        recorder = textwrap.dedent(recorder)
        (repo/'scripts/hpc_scheduler_record.py').write_text(recorder)
        (repo/'scripts/hpc_build_manifest.py').write_text(recorder)
        programs = {
            'module': '#!/bin/bash\nexit 0\n',
            'srun': '''
    import os,sys
    from pathlib import Path
    Path(os.environ['FLAG']).write_text('0')
    if sys.argv[-1]=='hostname': print('node-a\\nnode-b')
    ''',
            'mpiexec': '''
    import os,signal,time
    from pathlib import Path
    Path(os.environ['FLAG']).write_text('0')
    signal.signal(signal.SIGUSR1, lambda *args: None)
    Path(os.environ['IGA_CHECKPOINT_ROOT'],'manifest').write_text('test accepted checkpoint')
    time.sleep(.1);os.kill(os.getppid(),signal.SIGUSR1);time.sleep(.2)
    ''',
            'scontrol': '''
    import os,sys
    from pathlib import Path
    flag=Path(os.environ['FLAG'])
    with Path(os.environ['CONTROL_LOG']).open('a') as out: out.write(' '.join(sys.argv[1:])+'\\n')
    if sys.argv[1]=='update': flag.write_text('1')
    elif flag.read_text()!='1': print('requeue disabled',file=sys.stderr);sys.exit(2)
    '''
        }
        for name, code in programs.items():
            code = textwrap.dedent(code)
            path=binaries/name
            path.write_text(code if code.startswith('#!') else '#!'+sys.executable+'\n'+code)
            path.chmod(0o755)
        environment={k:v for k,v in os.environ.items() if not k.startswith(('SLURM_','IGA_','BASH_FUNC_','BASH_ENV'))}
        environment.update(PATH=str(binaries)+os.pathsep+os.environ['PATH'], FLAG=str(directory/'flag'),
            CONTROL_LOG=str(directory/'control.log'), IGA_REPO_ROOT=str(repo),
            IGA_GRAPH_CASE=str(directory/'case'), IGA_OUTPUT_ROOT=str(directory/'output'),
            IGA_CHECKPOINT_ROOT=str(directory/'checkpoint'), PETSC_DIR=str(directory/'petsc'),
            SLURM_JOB_ID='123', SLURM_JOB_NUM_NODES='2', SLURM_NTASKS='2',
            SLURM_NTASKS_PER_NODE='1', SLURM_CPUS_PER_TASK='1', LOCAL=str(directory/'local'),
            OMP_NUM_THREADS='1')
        with (directory/'stdout').open('w') as out, (directory/'stderr').open('w') as err:
            result=subprocess.run(['bash',str(directory/'spooled-script')],env=environment,stdout=out,stderr=err,timeout=20)
        record=json.loads((directory/'output/attempt-0/scheduler.json').read_text())
        results[variant]={'returncode':result.returncode,'scheduler':record,'control':(directory/'control.log').read_text()}
    return results['actual']


class GraphRequeueTests(unittest.TestCase):
    def test_launcher_flag_reset_is_corrected_after_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            result = run_wrapper(Path(directory))
        self.assertEqual(result['returncode'], 0, result)
        self.assertEqual(result['scheduler']['status'], 'checkpointed', result)
        self.assertTrue(result['scheduler']['checkpoint_requested'], result)
        self.assertTrue(result['control'].endswith('update JobId=123 Requeue=1\nrequeue 123\n'), result)


if __name__ == '__main__':
    unittest.main()
