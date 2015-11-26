import subprocess
import re
import matplotlib.pyplot as plt
import matplotlib.ticker

if __name__=='__main__':
    proc = subprocess.Popen('./fast_hadamard_transform', stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    xs, ys = [], []
    for line in stdout.split('\n'):
        mo = re.match('\[.*\(sz=(\d+)\)\].*? ([0-9.]+) us/iter.', line)
        if mo:
            sz = int(mo.group(1))
            us = float(mo.group(2))
            xs.append(sz)
            ys.append(us)

    fig, ax = plt.subplots(1, 1)
    ax.plot(xs, ys, 'o-')
    ax.set_title('Fast Hadamard Transform in C++')
    ax.set_xlabel('size')
    ax.set_ylabel('us/iter')
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
    plt.show()

