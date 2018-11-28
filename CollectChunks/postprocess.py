import os 
import sys
directory = sys.argv[1]

files = os.listdir(directory)
no_sucesses = 0
no_failures = 0
for f in files :
	if '.txt' in f :
		with open("./out/"+f) as name :
			content = name.readlines()
		if  len(content) != 0 :
			if 'successful' in content[0] : 
				no_sucesses += 1
			else :
				print 'failure in {0}'.format(f)
				no_failures += 1
		else :
			no_failures += 1
			print 'failure in {0}'.format(f)
print '-'*100
		
print "No of successes {0}".format(no_sucesses)
print "No of failures {0}".format(no_failures)
