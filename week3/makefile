v.png: plotter.py output.txt
	python plotter.py
#	#rm -f output.txt placas

placas: placas_C.c
	mpicc -o placas placas_C.c -lm
