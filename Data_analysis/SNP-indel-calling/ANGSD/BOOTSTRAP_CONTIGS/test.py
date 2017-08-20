from time import sleep
from sys import stderr

for i in range(10):
    print >> stderr, str(i)
    sleep(10)

