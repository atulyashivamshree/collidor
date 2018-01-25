import sys
import os
from cStringIO import StringIO
import argparse

RED = "\033[1;31m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"

# skips these items in the check
FILTERS = '-runtime/references,-build/include_what_you_use, -runtime/threadsafe_fn, -readability/casting, -build/include'
# Extensions to run cpplint on
EXTENSIONS = ['.h', '.cpp', '.cu']
# Folders to run cpplint recursively on
FOLDERS_TO_CHECK = ['src']


class FileList:

    def __init__(self, locations):
        self.files = []
        for loc in locations:
            os.path.walk(loc, self.step, '')

    def step(self, e, dirname, names):
        for name in names:
            for ext in EXTENSIONS:
                if name.endswith(ext):
                    self.files.append(dirname + '/' + name)

# returns the number of errors found and their description as a list


def testFile(file):
    result = [0, [], []]

    cpplint._cpplint_state.ResetErrorCounts()
    cpplint._SetFilters(FILTERS)

    original_stdout = sys.stdout
    original_stderr = sys.stderr
    temp_stdout = StringIO()
    temp_stderr = StringIO()
    sys.stdout = temp_stdout
    sys.stderr = temp_stderr
    cpplint.ProcessFile(file, 1)
    sys.stdout = original_stdout
    sys.stderr = original_stderr

    lines = temp_stdout.getvalue().split('\n')
    lines_err = temp_stderr.getvalue().split('\n')
    result[0] = cpplint._cpplint_state.error_count
    result[1] = lines
    result[2] = lines_err

    return result


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Runs cpplint on a set of folders recursively')

    parser.add_argument('loc', metavar='STYLEGUIDE_LOCATION', nargs=1,
                        help='relative location of the google styleguide folder')
    parser.add_argument('--batch', action='store_true',
                        help='Use this flag to run the program without pause')
    args = parser.parse_args()

    loc = args.loc[0]
    if not loc.endswith('/'):
        loc += '/'

    sys.path.append(loc + 'styleguide/cpplint')

    try:
        import cpplint
    except ImportError as err:
        print("Could not find the cpplint module in %sstyleguide/cpplint"%loc)
        print err
        exit()


    tree = FileList(FOLDERS_TO_CHECK)

    overall_results = []
    error_count = 0

    print "Run on source files under these folders : "
    print FOLDERS_TO_CHECK

    for file in tree.files:
        res = testFile(file)
        overall_results.append((len(overall_results), res[0]))
        error_count += res[0]

        if res[0] > 0:
            print '-------------'
            for err_msg in res[2]:
                print err_msg
            sys.stdout.write(RED)
            print 'Checking %s, Errors Found: %d' % (file, res[0])
            sys.stdout.write(RESET)
            if args.batch:
                continue
            ch = raw_input('Enter s to skip, antyhing else to quit\n')
            if ch == 's':
                continue
            else:
                break
        sys.stdout.write(GREEN)
        print 'Checking %s, Errors Found: %d' % (file, res[0])
        sys.stdout.write(RESET)

    print '\n\n==================\n      Summary\n'
    print 'Total files checked %d, Errors found: %d' %\
        (len(overall_results), error_count)
