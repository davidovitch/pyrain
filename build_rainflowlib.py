"""
f2py c_library.pyf c_functions.c -c

See also http://www.scipy.org/Cookbook/CompilingExtensionsOnWindowsWithMinGW
"""
import os

def which(program):
    """
    Test if program exists
    ======================
    
    In order to test if a certain executable exists, it will search for the 
    program name in the environment variables.
    If program is a full path to an executable, it will check it exists
    
    Copied from:
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/
    It is supposed to mimic the UNIX command "which"
    """
    
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def compile_all():
    
    # regardless of platform, try to figure out which f2py call is in the path
    # define possible options
    f2py_call_list = ('f2py','f2py2.6','f2py2.7','f2py.py', 'f2py2.5')
    
    no_f2py = True
    for k in f2py_call_list:
        # if the call command exists in the path, it will return the path as
        # a string, otherwise it will return None
        f2py_path = which(k)
        if not f2py_path:
            # didn't find the current call k, continue looking
            pass
        else:
            # current call k is in the path
            f2py_call = k
            no_f2py = False
            break
    
    # raise exception if f2py is not found
    if no_f2py:
        raise UserWarning, \
        'Couldn\'t locate f2py. Should be part of NumPy installation.'
    else:
        print '='*75
        print 'compiling rainflowlib.c'
        print '='*75
        print 'found f2py in:', f2py_path
    
    # on Windows: Install microsoft visual c++ .NET 2003 to run the following 
    # build command
    # on posix: install gcc and gfortran
    compile_format = f2py_call + ' %s %s -c'
    
    pyfs = ('rainflowlib.pyf',)
    files =('rainflowlib.c',)

    for pyf,file in zip(pyfs,files):
        os.system(compile_format % (pyf,file))

if __name__=='__main__':
    compile_all()
