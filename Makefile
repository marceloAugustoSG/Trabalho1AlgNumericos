all: programa

programa: main.c funcoes.c funcoes.h
	@echo "Compilando o programa"
	gcc -o main main.c funcoes.c -Wall -Wextra -Wno-unused-result -Wpedantic -O0

clean:
	rm -rf *.o *.gch

teste: programa
	@echo "Executando teste..."
	./programa < input/ex1_3x3.in