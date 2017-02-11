PNG = Potential_?.png

all: ${PNG} placas

placas: placas_C.c
	mpicc -o placas placas_C.c -lm

${PNG}: plotter.py output.txt
#	mv camilor_placas.o????? output.txt
	python plotter.py
#	#rm -f output.txt placas

#output.txt: camilor_placas.o?????
#	mv camilor_placas.o????? output.txt

#camilor_placas.o?????: submit_job.sh
#	qsub submit_job.sh
