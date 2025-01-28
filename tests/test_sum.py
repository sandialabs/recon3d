import python_template.sum as sm


def test_sum():
    """Tests the sum function"""
    a = 2
    b = 3
    known_val = 5
    found_val = sm.calc_sum(a=a, b=b)
    assert found_val == known_val
