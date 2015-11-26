#!env python
# Fast Hadamard Transform
import numpy as np

def _power_of_two(sz):
    n = 0
    while 2**n < sz:
        n += 1
    return n

def fhtpad(x):
    u'''pad 1-d array x to the minimal power-of-two'''
    x = x.ravel()
    n = _power_of_two(x.size)
    if x.size < 2**n:
        tmp = np.zeros(2**n)
        tmp[:x.size] = x
        x = tmp
    return x

def fht(x, unitary=False):
    u'''Fast Hadamard Transform of 2^n size array x.'''
    x = x.ravel()
    n = _power_of_two(x.size)
    assert x.size == 2**n
    t0, t1 = np.array(x), np.zeros_like(x)
    for step in [2**(n - i - 1) for i in range(n)]:
        skip = step*2
        for j in range(0, t0.size, skip):
            t1[j       :j + step] = t0[j:j + step] + t0[j + step:j + skip]
            t1[j + step:j + skip] = t0[j:j + step] - t0[j + step:j + skip]
        t0, t1 = t1, t0
    if unitary:
        return t0/(np.sqrt(2.0)**n)
    else:
        return t0

def ifht(x, unitary=False):
    u'''Inverse Fast Hadamard Transform of 2^n size array x.'''
    if unitary:
        # in unitary mode, fht(x) == ift(x)
        return fht(x, True)
    else:
        # in non-unitary mode, ifht(fht(x)) == x but not fht(x) != ift(x).
        n = _power_of_two(x.size)
        return fht(x)/2.0**n

if __name__=='__main__':

    # non-unitary
    arr = np.array([1, 0, 1, 0, 0, 1, 1], np.float32) # size=7
    arr_ht = fht(fhtpad(arr))
    arr_ht_ht = fht(arr_ht)
    arr_ht_iht = ifht(arr_ht)
    print 'Non-unitary..'
    print 'x         :', arr
    print 'HT(x)     :', arr_ht
    print 'HT(HT(x)) :', arr_ht_ht
    print 'IHT(HT(x)):', arr_ht_iht
    assert np.allclose(fhtpad(arr), arr_ht_iht)

    # unitary
    arr = np.array([1, 0, 1, 0, 0, 1, 1, 0], np.float32) # size=8
    arr_ht = fht(arr, True)
    arr_ht_ht = fht(arr_ht, True)
    arr_ht_iht = ifht(arr_ht, True)
    print 'Unitary..'
    print 'x         :', arr
    print 'HT(x)     :', arr_ht
    print 'HT(HT(x)) :', arr_ht_ht
    print 'IHT(HT(x)):', arr_ht_iht
    assert np.allclose(arr, arr_ht_iht)

