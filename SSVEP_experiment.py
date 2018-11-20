import random
import time

import SerialPort as Sp


def task_1(sp):
	freqs = dict.fromkeys([4, 7, 10, 13, 16, 20, 25, 30, 35, 40], 0)

	times = [3, 4, 5, 6, 7]

	time_length = 5

	prev_freq = 0

	number_iterations = 1

	while freqs:
		print '[LOG] frequencies left', len(freqs.values())

		_t = random.choice(times)
		print '[LOG] sleeping for', _t
		time.sleep(_t)

		freq = random.choice(freqs.keys())

		while freq == prev_freq:
			freq = random.choice(freqs.keys())

		freqs[freq] += 1

		print '[LOG] current frequency', freq

		sp.blinkSSVEP([freq, 0], 1, 1)

		time.sleep(time_length)

		sp.blinkSSVEP([0, 0], 1, 1)

		if freqs[freq] >= number_iterations:
			del freqs[freq]
			print '[LOG] removing frequency', freq


sp = Sp.SerialPort('/dev/ttyUSB0')
sp.open()

sp.blinkSSVEP([0, 0], 1, 1)

print('Starting Task 1')
raw_input('Are you ready? ...')
task_1(sp)
raw_input('Are you ready? ...')
task_1(sp)

sp.close()
