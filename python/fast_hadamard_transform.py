#!env python
# Fast Hadamard Transform
import numpy as np
import scipy.misc
import scipy.ndimage
import matplotlib
import matplotlib.pyplot as plt

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

def fht2(x, unitary=False):
    u'''2D Fast Hadamard Transform'''
    x = np.array([fht(row, unitary) for row in x])
    x = np.array([fht(col, unitary) for col in x.T])
    return x

def ifht2(x, unitary=False):
    u'''2D Inverse Fast Hadamard Transform'''
    x = np.array([ifht(row, unitary) for row in x])
    x = np.array([ifht(col, unitary) for col in x.T])
    return x

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

    # 2D
    img = scipy.misc.lena()
    img = scipy.ndimage.zoom(img, 1.0/8)
    matplotlib.rc('font', size=9)
    fig, axs = plt.subplots(4, 3, figsize=(12, 9))

    img_ht = fht2(img)
    img_ht_iht = ifht2(img_ht)
    assert np.allclose(img, img_ht_iht)

    ax = axs[0, 0]; ax.imshow(img, cmap='gray'); ax.set_title('org')
    ax = axs[0, 1]; ax.imshow(img_ht); ax.set_title('HT(img)')
    ax = axs[0, 2]; ax.imshow(img_ht_iht, cmap='gray'); ax.set_title('IHT(HT(img))')
    ax = axs[1, 0]; ax.hist(img.ravel(), bins=50, edgecolor='none'); ax.set_title('org')
    ax = axs[1, 1]; ax.hist(img_ht.ravel(), bins=50, edgecolor='none'); ax.set_title('HT(img)')
    fig.delaxes(axs[1, 2])

    img_u_ht = fht2(img, True)
    img_u_ht_iht = ifht2(img_u_ht, True)
    assert np.allclose(img, img_u_ht_iht)

    ax = axs[2, 0]; ax.imshow(img, cmap='gray'); ax.set_title('org')
    ax = axs[2, 1]; ax.imshow(img_u_ht); ax.set_title('uHT(img)')
    ax = axs[2, 2]; ax.imshow(img_u_ht_iht, cmap='gray'); ax.set_title('uIHT(uHT(img))')
    ax = axs[3, 0]; ax.hist(img.ravel(), bins=50, edgecolor='none'); ax.set_title('org')
    ax = axs[3, 1]; ax.hist(img_u_ht.ravel(), bins=50, edgecolor='none'); ax.set_title('uHT(img)')
    fig.delaxes(axs[3, 2])

    fig.subplots_adjust(hspace=0.4)
    fig.suptitle('Fast Hadamard Transform of 2D image', fontsize=12)
    plt.show()

