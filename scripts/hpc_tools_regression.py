#!/usr/bin/env python3
"""Validate auxiliary MPI tools against saved pre-change binaries and controlled faults."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import time
from hpc_inventory import digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-source',type=Path,required=True)
    parser.add_argument('--baseline-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--launcher',default='mpiexec --oversubscribe')
    args=parser.parse_args(); repo=Path(__file__).resolve().parents[1]
    out=args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=False)
    fixture=out/'fixture'; shutil.copytree(args.fixture_source,fixture)
    # This bounded one-element v5 test asset uses three outlet labels; the
    # legacy conversion contract only defines 0/1/2. Adapt the test copy only.
    for filename in ('serial.ntiga','group.ntiga'):
        path=fixture/filename; data=bytearray(path.read_bytes())
        assert struct.unpack_from('=I',data,8)[0]==5
        assert struct.unpack_from('=QQ',data,16)==(1,64)
        offset=struct.unpack_from('=Q',data,88)[0]
        assert struct.unpack_from('=I',data,offset+16)[0]==64
        labels=list(struct.unpack_from('=6i',data,offset+20))
        assert 3 in labels
        struct.pack_into('=6i',data,offset+20,*[2 if label==3 else label for label in labels])
        path.write_bytes(data)
    mesh=fixture/'controlmesh.vtk'; before,labels=mesh.read_text().rsplit('LOOKUP_TABLE default',1)
    labels=labels.split(); assert len(labels)==64
    mesh.write_text(before+'LOOKUP_TABLE default\n'+'\n'.join('2' if x=='3' else x for x in labels)+'\n')
    tools=('iga_mesh_check','iga_assembly_smoke','iga_transport')
    current={n:repo/'solvers/cpu'/n for n in tools}
    old={n:args.baseline_dir.resolve()/n for n in tools}
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
             PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps')
    summary=dict(status='running',binaries={str(p):digest(p) for p in [*current.values(),*old.values()]},
                 fixture={str(p):digest(p) for p in fixture.iterdir() if p.is_file()},cases=[],comparisons=[])
    def run(name,tool,children,expected=0,stage='',baseline=False):
        directory=out/name; directory.mkdir()
        binary=(old if baseline else current)[tool]
        cmd=['timeout','--kill-after=5s','90s',*shlex.split(args.launcher)]
        for i,(settings,argv) in enumerate(children):
            if i: cmd += [':']
            cmd += ['-np','1','env',*[k+'='+v for k,v in settings.items()],str(binary),*map(str,argv)]
        record=dict(name=name,tool=tool,baseline=baseline,command_argv=cmd,expected=expected,stage=stage,status='running')
        summary['cases'].append(record); started=time.monotonic(); log=directory/'launcher.log'
        with log.open('x') as stream: result=subprocess.run(cmd,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT)
        text=log.read_text(); record.update(returncode=result.returncode,elapsed_s=time.monotonic()-started,log_sha256=digest(log))
        if result.returncode!=expected or stage not in text: raise RuntimeError(f'{name}: expected exit {expected} and {stage!r}')
        if expected==1 and any(x in text for x in ('transport_v2 nodes=','global_rows=','minimum_detJ=')):
            raise RuntimeError(f'{name}: rejected input published success summary')
        record['status']='passed';print(name,'passed',flush=True)
        return text
    def arguments(tool,db,ranks=2,target=None):
        if tool=='iga_mesh_check': return [db]
        if tool=='iga_assembly_smoke': return [db,'2']
        return [db,fixture,'2',target or out/f'field-{ranks}.txt']
    try:
        for ranks,db in ((1,fixture/'serial.ntiga'),(2,fixture/'group.ntiga')):
            for tool in tools:
                texts=[]
                for baseline in (True,False):
                    tag=('before' if baseline else 'after')+f'-{tool}-{ranks}'
                    target=out/(tag+'.txt'); argv=arguments(tool,db,ranks,target)
                    texts.append(run(tag,tool,[({},argv)]*ranks,baseline=baseline))
                if tool!='iga_transport':
                    prefix='elements=' if tool=='iga_mesh_check' else 'global_rows='
                    a=[x for x in texts[0].splitlines() if x.startswith(prefix)]
                    b=[x for x in texts[1].splitlines() if x.startswith(prefix)]
                    if a!=b or len(a)!=1: raise RuntimeError('geometry/assembly summary changed')
                    summary['comparisons'].append(dict(tool=tool,ranks=ranks,summary=a[0],exact=True))
                else:
                    a=out/f'before-{tool}-{ranks}.txt';b=out/f'after-{tool}-{ranks}.txt'
                    if a.read_bytes()!=b.read_bytes(): raise RuntimeError('legacy field changed within same rank count')
                    values=[float(x) for line in b.read_text().splitlines() for x in line.split()[1:]]
                    if len(values)!=128 or not all(math.isfinite(x) for x in values) or math.hypot(*values)==0:
                        raise RuntimeError('invalid or zero legacy reference')
                    summary['comparisons'].append(dict(tool=tool,ranks=ranks,field_sha256=digest(b),exact=True))
        fields=[]
        for ranks in (1,2):
            fields.append([float(x) for row in (out/f'after-iga_transport-{ranks}.txt').read_text().splitlines() for x in row.split()[1:]])
        error=math.hypot(*(a-b for a,b in zip(*fields)))/math.hypot(*fields[0])
        if error>1e-6: raise RuntimeError('legacy serial/MPI relative L2 failed')
        summary['serial_mpi_relative_l2']=error
        for tool in tools:
            stage={'iga_mesh_check':'mesh check database input','iga_assembly_smoke':'assembly smoke input','iga_transport':'legacy transport database input'}[tool]
            for mode in ('missing','partitions','arguments','threads'):
                children=[]
                for rank in range(2):
                    db=fixture/('serial.ntiga' if mode=='partitions' else 'group.ntiga')
                    if mode=='missing' and rank==1: db=fixture/'missing.ntiga'
                    argv=arguments(tool,db);settings={}
                    if mode=='arguments' and rank==1: argv=[]
                    if mode=='threads' and rank==1: settings={'OMP_NUM_THREADS':'bad'}
                    children.append((settings,argv))
                run(tool+'-'+mode,tool,children,1,'execution resource preflight' if mode=='threads' else stage)
        for value in ('0','-1','2junk','4294967297','9223372036854775808'):
            run('assembly-fields-'+value,'iga_assembly_smoke',[({},[fixture/'group.ntiga','2']),({},[fixture/'group.ntiga',value])],1,'assembly smoke input')
        run('assembly-fields-different','iga_assembly_smoke',[({},[fixture/'group.ntiga',n]) for n in ('2','3')],1,'assembly smoke field agreement')
        for mode in ('missing-case','bad-steps','different-steps','different-output','output-directory','output-full','preonly'):
            children=[]
            for rank in range(2):
                argv=arguments('iga_transport',fixture/'group.ntiga');settings={}
                if mode=='missing-case' and rank==1: argv[1]=fixture/'missing'
                if mode=='bad-steps' and rank==1: argv[2]='2x'
                if mode=='different-steps' and rank==1: argv[2]='1'
                if mode=='different-output' and rank==1: argv=argv[:3]
                if mode=='output-directory': argv[3]=fixture
                if mode=='output-full': argv[3]='/dev/full'
                if mode=='preonly': settings={'PETSC_OPTIONS':'-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'}
                children.append((settings,argv))
            stage='asset content read' if mode=='missing-case' else 'legacy transport case input' if mode=='bad-steps' else 'legacy transport execution agreement' if mode.startswith('different') else 'legacy transport linear solve' if mode=='preonly' else 'legacy transport output'
            run('legacy-'+mode,'iga_transport',children,1,stage)
        flat=out/'flat.ntiga';data=bytearray((fixture/'group.ntiga').read_bytes())
        end=struct.unpack_from('=Q',data,40)[0]
        for i in range(64): struct.pack_into('=d',data,end-64*24+i*24+16,0.0)
        flat.write_bytes(data)
        for baseline in (True,False):
            text=run('flat-'+('before' if baseline else 'after'),'iga_mesh_check',[({},[flat])]*2,2,'bad_elements=1 bad_samples=64',baseline)
            if 'minimum_detJ=0' not in text: raise RuntimeError('bad Jacobian semantics changed')
        summary['status']='passed'
    except BaseException:
        summary['status']='failed';raise
    finally:
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__': main()
