"""
not a testsuite, actually, but a program that prints instructions
to manually test the GUI
"""

import sys

if __name__ == '__main__':

    sys.stdout.write('show help...');
    sys.stdin.readline()

    sys.stdout.write('scan single file...');
    sys.stdin.readline()
    sys.stdout.write('view results...');
    sys.stdin.readline()
    sys.stdout.write('store results...');
    sys.stdin.readline()

    sys.stdout.write('scan paired file...');
    sys.stdin.readline()
    sys.stdout.write('view results...');
    sys.stdin.readline()

    sys.stdout.write('scan multiple files...');
    sys.stdin.readline()
    sys.stdout.write('view results...');
    sys.stdin.readline()
    sys.stdout.write('store results...');
    sys.stdin.readline()

    sys.stdout.write('explore results...');
    sys.stdin.readline()
    sys.stdout.write('summarize...');
    sys.stdin.readline()

    sys.stdout.write('scan broken file...');
    sys.stdin.readline()

    sys.stdout.write('show log...');
    sys.stdin.readline()

    sys.stdout.write('done: close window\n\n')

