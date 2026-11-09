"""Verify a wheel in a fresh environment outside the source checkout.

The core-only job deliberately omits legacy NetworkX metadata dependencies;
ordinary installations continue to include them for backward compatibility.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import venv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wheel', required=True)
    parser.add_argument('--constraints', required=True)
    parser.add_argument('--core-only', action='store_true')
    parser.add_argument('--report', required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    wheel, constraints = Path(args.wheel).resolve(), Path(args.constraints).resolve()
    report = Path(args.report).resolve()
    with tempfile.TemporaryDirectory(prefix='pmg-wheel-') as name:
        outside = Path(name)
        venv.EnvBuilder(with_pip=True).create(outside/'venv')
        python = outside/'venv'/('Scripts/python.exe' if sys.platform == 'win32' else 'bin/python')

        def run(*command, capture=False):
            result = subprocess.run([str(python), *command], cwd=outside,
                                    text=True, capture_output=capture, check=True)
            return result.stdout if capture else None

        dependencies = ['numpy','sympy','pytest']
        if not args.core_only:
            dependencies += ['networkx','pandas']
        run('-m','pip','install','-c',str(constraints),*dependencies)
        run('-m','pip','install','--no-deps',str(wheel))
        for folder in ('tests','examples'):
            shutil.copytree(root/folder, outside/folder, ignore=shutil.ignore_patterns('__pycache__'))
        check = ('import json,sys,polyform_algebra as p; '
                 'from pathlib import Path; '
                 'assert Path(p.__file__).is_relative_to(Path(sys.prefix)); '
                 'assert p.get_version()==p.__version__; '
                 'print(json.dumps(dict(version=p.get_version(),python=sys.version,origin=p.__file__)))')
        metadata = json.loads(run('-c',check,capture=True))
        if args.core_only:
            run('-c', "import importlib.util; assert all(importlib.util.find_spec(x) is None for x in ('networkx','pandas','scipy'))")
            targets = ['tests/test_core_data.py','tests/test_core_network_math.py',
                       'tests/test_core_reductions.py','tests/test_core_algebra.py',
                       'tests/test_delivery.py','tests/test_connector_contracts.py','tests/test_comment_files.py','tests/reference']
        else:
            run('-m','pip','check')
            targets = ['tests']
        try:
            run('-m','pytest','-q',*targets,'--junitxml=result.xml')
        finally:
            # Preserve evidence for failing jobs as well as successful ones.
            report.parent.mkdir(parents=True,exist_ok=True)
            metadata.update(core_only=args.core_only, freeze=run('-m','pip','freeze',capture=True))
            report.write_text(json.dumps(metadata,indent=2),encoding='utf-8')
            if (outside/'result.xml').exists():
                shutil.copyfile(outside/'result.xml',report.with_suffix('.xml'))


if __name__ == '__main__':
    main()
