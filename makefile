cc = g++
prom = calc
deps = Parameters.h Solver.h
obj = main.o Parameters.o Solver.o

$(prom): $(obj)
		   $(cc) -o $(prom) $(obj)

%.o: %.cpp $(deps)
	$(cc) -c $< -o $@

clean:
	rm -rf $(prom) $(obj) *.dat
