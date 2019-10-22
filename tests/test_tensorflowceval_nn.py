import os.path
import pytest
from bamboo.root import findLibrary
pytest.mark.skipif(not findLibrary("libtensorflow"))

testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def tfc_nn1():
    from bamboo.root import gbl, loadTensorflowC
    loadTensorflowC()
    yield gbl.bamboo.TensorflowCEvaluator(os.path.join(testData, "nn1.pb"), "dense_1_input", "dense_3/Softmax", 5, 2)

def test_tfc_load(tfc_nn1):
    assert tfc_nn1

def test_tfc_evalZero(tfc_nn1):
    from bamboo.root import gbl
    v_in = getattr(gbl, "std::vector<float>")()
    for i in range(5):
        v_in.push_back(.7)
    v_out = tfc_nn1(v_in)
    assert v_out.size() == 2

def test_tfc_evalRandom(tfc_nn1):
    from bamboo.root import gbl
    import random
    v_in = getattr(gbl, "std::vector<float>")()
    for i in range(5):
        v_in.push_back(random.uniform(0., 1.))
    v_out = tfc_nn1(v_in)
    assert v_out.size() == 2
