__author__ = 'Stephen G. Gaffney'


dbvars = {'host': 'localhost', 'db': 'CancerDB',
          'read_default_file': "~/.my.cnf"}


class NoIntervalsException(Exception):
    pass


class LookupFailedException(Exception):
    pass



