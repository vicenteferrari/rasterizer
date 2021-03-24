CC=g++

main_source := $(wildcard src/*.cpp)
main_headers := $(wildcard src/*.h)
data := $(wildcard data/*)

.PHONY: all
all: rasterizer.exe

rasterizer.exe: rasterizer.o
	$(CC) rasterizer.o -o rasterizer.exe -lobjloader -lmingw32 -lSDL2main -lm -lpthread -lSDL2

rasterizer.o: $(main_source) $(main_headers)
	$(CC) -g -c src/rasterizer.cpp -o rasterizer.o


.PHONY: clean
clean:
	rm -f rasterizer.o
	rm -f rasterizer.exe
