"""Exact symbolic identifiability witness; no trajectory simulation."""
import json,sys
from pathlib import Path
import sympy as sp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.support_information import POWERS


def main():
    x,v=sp.symbols('x v');rhs=-x+v-x*x*v
    functions=[x**int(a)*v**int(b) for a,b in POWERS];rows=[]
    for _ in range(9):
        rows.append([f.subs({x:2,v:0}) for f in functions])
        functions=[sp.expand(v*sp.diff(f,x)+rhs*sp.diff(f,v)) for f in functions]
    matrix=sp.Matrix(rows);det=matrix.det();assert det==-65265167726360395776
    a,b,c,d,delta=sp.symbols('a b c d delta')
    energy=v*v/2-b*x*x/(2*a)-delta*x**4/(4*a)
    derivative=sp.expand(sp.diff(energy,x)*a*v+sp.diff(energy,v)*(b*x+c*v+d*x*x*v+delta*x**3))
    assert sp.simplify(derivative-(c+d*x*x)*v*v)==0
    report={'arithmetic':'exact symbolic integers/rationals','powers':POWERS.tolist(),
        'jet_rows':[[int(q) for q in row] for row in rows],'determinant':str(det),'rank':int(matrix.rank()),
        'energy_derivative':str(derivative),'trajectory_simulations':False}
    (ROOT/'reports/vdp_exact_identifiability.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['determinant','rank','energy_derivative']}))


if __name__=='__main__':main()
