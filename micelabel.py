# the following doesn't work, because python gets the columnnames wrong
#
#python -c "
#import sys
#import csv
#
#input=\"$1\"
#
#checkLabelNr=input.isdigit()
#
#f = open(\"$DEFS\", 'rt')
#try:
#    reader = csv.DictReader(f)
#    colnames = reader.fieldnames
#    for row in reader:
#		if checkLabelNr:
#			if input in row.values():
#				print row
#finally:
#    f.close()
#
#
#print colnames
#sys.exit(1)
#"
#

