from subprocess import Popen, PIPE
from os import listdir, kill
from time import sleep

exeFile = '/home/users/bybeep/Coursera/check_brackets_in_code/dist/Debug/GNU-Linux/check_brackets_in_code'

testDIR = '/home/users/bybeep/Coursera/check_brackets_in_code/tests'

testFiles = listdir(testDIR)


for fd in testFiles:
    if fd[-1] == 'a':
        continue
    p = Popen(exeFile, shell=True, stdout=PIPE, stdin=PIPE)
    # sleep(0.1)
    print("testcase: " + fd,)
    with open(testDIR+'/'+fd) as inputData:
        content = inputData.read()

    with open(testDIR+'/'+fd + '.a') as answer:
        answerStr = answer.read()

    p.stdin.write(content)
    # sleep(0.1)
    p.stdin.flush()
    result = p.stdout.readlines()
    # sleep(0.1)
    p.stdout.flush()
    # print(result)

    kill(p.pid,9)

    sleep(0.1)

    if (result[0] == answerStr.strip()):
        print(" PASSED")

