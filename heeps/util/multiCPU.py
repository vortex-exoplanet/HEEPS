import numpy as np
import multiprocessing as mpro
from functools import partial
from sys import platform
import time

def multiCPU(f, posargs=[], posvars=[], kwargs={}, multi_out=False,
        estimate_time=False, case=None, cpu_count=None, verbose=True):

    assert type(posargs)==list, 'posargs must be a list [] of positionnal arguments'
    assert type(posvars)==list, 'posvars must be a list [] of positionnal variables for multiprocessing'
    assert type(kwargs)==dict, 'kwargs must be a dict {} of keyword arguments'

    # format case
    case = '' if case is None else case + ' '
    # get initial time
    t0 = time.time()

    # run one and time it
    if estimate_time is True:
        x0 = list(x[0] for x in posvars)
        args = posargs + x0
        res_one = f(*args, **kwargs)
        t1 = time.time() - t0
        nsim = len(posvars[0])

    # run multicore
    if cpu_count != 1 and platform in ['linux', 'linux2', 'darwin']:
        if cpu_count == None:
            cpu_count = mpro.cpu_count()
        if verbose is True:
            print('   %s, %susing %s cores'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), case, cpu_count))
            if estimate_time is True:
                texp = t1*np.ceil(nsim/cpu_count)
                print('   expected time to complete: %s seconds'%np.round(texp, 2))
        p = mpro.Pool(cpu_count)
        func = partial(f, *posargs, **kwargs)
        res_multi = np.array(p.starmap(func, zip(*posvars)))
        if multi_out is True:
            res_multi = tuple(np.array(x) for x in res_multi.swapaxes(0, 1))
        p.close()
        p.join()

    # run singlecore
    else:
        if verbose is True:
            print('   %s, %susing 1 core'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), case))
            if estimate_time is True:
                texp = t1*nsim
                print('   expected time to complete: %s seconds'%np.round(texp, 2))
        
        for i, xi in enumerate(zip(*posvars)):
            args = posargs + list(xi)
            res = f(*args, **kwargs)
            if i == 0:
                res_multi = res
            else:
                if multi_out is True:
                    res_multi = tuple(anystack(x, y) for x, y in zip(res_multi, res))
                else:
                    res_multi = anystack(res_multi, res)
    
    if verbose is True:
        print('   %s, completed in %s seconds'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), np.round(time.time() - t0, 2)))

    return res_multi

def anystack(x,y):
    dsum = x.ndim + y.ndim
    if dsum < 2:
        z = np.hstack((x,y))
    elif dsum < 4:
        z = np.vstack((x,y))
    else:
        z = np.dstack((x.T,y.T)).T
    return z