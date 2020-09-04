def round2odd(x):

    '''
    x: float
        real number
    '''

    x1 = x + 0.5
    n = int(x1)
    if not (n % 2):
        n += 1 if int(x1) == int(x) else - 1
    
    return n