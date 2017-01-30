all: precision_float.c graph_precisiones.py
	gcc precision_float.c -o pre.x
	./pre.x
	python graph_precisiones.py
	
clean: 
	rm pre.x Precisiones_float.txt Precisiones_float.pdf
