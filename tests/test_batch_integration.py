import numpy as np
from scripts.run_revision_fast import batch_vdp
from src.systems import trajectory


def test_batched_vdp_matches_individual_integration():
    t=np.linspace(0,20,101);initials=np.array([[2.,0.],[1.8,.2],[2.2,-.1]])
    batch=batch_vdp(t,initials)
    reference=np.stack([trajectory('vdp',t,initial) for initial in initials],axis=1)
    np.testing.assert_allclose(batch,reference,rtol=1e-7,atol=1e-8)
