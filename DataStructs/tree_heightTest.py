from subprocess import Popen, PIPE
from os import listdir, kill
from time import sleep

exeFile = '/home/users/bybeep/Coursera/tree_height/dist/Debug/GNU-Linux/tree_height'

testDIR = '/home/users/bybeep/Coursera/tree_height/tests'

testFiles = listdir(testDIR)


for fd in testFiles:
    if fd[-1] == 'a':
        continue
    p = Popen(exeFile, shell=True, stdout=PIPE, stdin=PIPE)
    # sleep(0.1)
    print("testcase: " + fd)
    with open(testDIR+'/'+fd) as inputData:
        content = inputData.read()

    with open(testDIR+'/'+fd + '.a') as answer:
        answerStr = answer.read().replace("\r\n", " ")

    p.stdin.write(content)
    # sleep(0.1)
    p.stdin.flush()
    result = p.stdout.read().replace("\n", " ")
    # sleep(0.1)
    p.stdout.flush()
    # print(result)

    kill(p.pid,9)

    sleep(0.1)

    if (result == answerStr):
        print(" PASSED")

