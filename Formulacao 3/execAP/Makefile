CFLAGS= -Wall -m64 -g -w
#CFLAGS= -Wall -m64 -o2 
CXX=g++

ILOG= /opt/ibm/ILOG/CPLEX_Studio2211
CPPFLAGS= -DIL_STD -I$(ILOG)/cplex/include -I$(ILOG)/concert/include
CPLEXLIB=-L$(ILOG)/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L$(ILOG)/concert/lib/x86-64_linux/static_pic -lconcert -lm -lpthread -ldl

comp-c: 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -O2 -o form3condAP formulacao3AP.cpp  $(CPLEXLIB) 

comp-form3: 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o form3  formulacao3.cpp  $(CPLEXLIB) 



run:comp-mono
	./mono FL-50-50	
clean:
	rm -f  *.out *.aux *.log *.nav *.snm *.out *.toc 
