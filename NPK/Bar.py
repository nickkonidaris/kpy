
import sys


'''
    Create a console bar indicating the amount of time that's passed
'''


def setup(toolbar_width=40):
    ''' Initialize console with toolbar_width spaces between [ ]

    Args:
        toolbar_width: Number of spaces between the brackets

    Returns:
        toolbar_width [int]'''

    toolbar_width = 40
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))

    return toolbar_width


def update(char='-'):
    '''Prints char to indicate an update on the toolbar

    Args:
        char: The character to print'''
    sys.stdout.write(char)
    sys.stdout.flush()

def done():
    '''Carriage return and flush the console'''
    sys.stdout.write("\n")
    sys.stdout.flush()

 
