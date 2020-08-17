import os.path
import pytest
from helpers import isclose_float
# check if lwtnn is found (the bridge library should be available then)
from bamboo.root  import findLibrary
lwtnn_found = findLibrary("liblwtnn")
if not lwtnn_found:
    if "VIRTUAL_ENV" in os.environ:
        lwtnnlib_guessed_path = os.path.join(os.environ["VIRTUAL_ENV"], "lib", "liblwtnn.so")
        if os.path.exists(lwtnnlib_guessed_path):
            lwtnn_found = lwtnnlib_guessed_path
pytestmark = pytest.mark.skipif(not lwtnn_found, reason="lwtnn library not found")

testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def lwtnn_nn1():
    from bamboo.root import loadlwtnn, gbl
    loadlwtnn()
    inputs = getattr(gbl, "std::vector<std::pair<std::string,std::string>>")()
    input_item = getattr(gbl, "std::pair<std::string,std::string>")
    for i in range(5):
        inputs.push_back(input_item("node_0", "variable_{0:d}".format(i)))
    outputs = getattr(gbl, "std::vector<std::string>")()
    for i in range(2):
        outputs.push_back("out_{0:d}".format(i))
    yield gbl.bamboo.LwtnnEvaluator(os.path.join(testData, "nn1_lwtnn.json"), inputs, outputs)

def test_lwtnn_load(lwtnn_nn1):
    assert lwtnn_nn1

def test_lwtnn_evalZero(lwtnn_nn1):
    from bamboo.root import gbl
    v_in = getattr(gbl, "std::vector<float>")()
    for i in range(5):
        v_in.push_back(.7)
    v_out = lwtnn_nn1.evaluate(v_in)
    assert ( v_out.size() == 2 ) and isclose_float(sum(v_out), 1.)

def test_lwtnn_evalRandom(lwtnn_nn1):
    from bamboo.root import gbl
    import random
    v_in = getattr(gbl, "std::vector<float>")()
    for i in range(5):
        v_in.push_back(random.uniform(0., 1.))
    v_out = lwtnn_nn1.evaluate(v_in)
    assert ( v_out.size() == 2 ) and isclose_float(sum(v_out), 1.)
