'''
Module 'stats' - Contains functions for basic statistical operations
'''


###############
## FUNCTIONS ##
###############

def fraction(counts, total):
    '''
    Fraction calculator
    
    Input:
    1. counts: numerator
    2. total: denominator
    
    Output:
    1. fract: division result    
    '''
    if total > 0:
        fract = counts/total
    else:
        fract = None
    
    return fract


