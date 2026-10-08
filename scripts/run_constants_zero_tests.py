#!/usr/bin/env python3
from pathlib import Path
import argparse,subprocess,shutil,json
ap=argparse.ArgumentParser();ap.add_argument('--library-root',type=Path,required=True);ap.add_argument('--label',default='constants_zero_tests')
a=ap.parse_args();root=Path(__file__).resolve().parent;lib=a.library_root.resolve();out=lib/'ifx'/a.label;out.mkdir(exist_ok=True)
bin=out/'driver';r=subprocess.run(['mpiifx','-O2','-I'+str(lib/'ifx/mod'),str(lib/'test/driver_constants_zero.f90'),str(lib/'ifx/libSpeedCHEM64.a'),'-o',str(bin)],capture_output=True,text=True);(out/'compile.log').write_text(r.stdout+r.stderr);r.check_returncode()
summary={}
for units in ['CAL','KCAL','JOULES','KJOULES','KELV']:
 d=out/units;d.mkdir(exist_ok=True)
 shutil.copy2(lib/'test/data_plog/therm.dat',d/'therm.dat')
 (d/'chem.inp').write_text(f'''ELEMENTS
H O N
END
SPECIES
H2 O2 H2O H O OH HO2 H2O2 N2
END
REACTIONS {units}
H2 + O2 = OH + OH 1e10 0 10000
REV / 1e10 0 10000 /
H + O + M => OH + M 1e10 0 10000
N2/ 0 / H2O/ 2 /
H + O2 (+M) => HO2 (+M) 1e10 0 10000
LOW / 1e15 0 10000 /
O + H2 => OH + H 1e10 0 10000
PLOG / 0.1 1e10 0 10000 /
PLOG / 10 1e10 0 10000 /
END
''')
 r=subprocess.run([str(bin),units],cwd=d,capture_output=True,text=True,timeout=30);(d/'driver.log').write_text(r.stdout+r.stderr);summary[units]={'exit_code':r.returncode,'pass':r.returncode==0 and 'RESULT: PASS' in r.stdout}
print(json.dumps(summary,indent=2));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if not all(x['pass'] for x in summary.values()):raise SystemExit(2)
