
.SUFFIXES: .o .C .x

OBJ=create_kpp_module.o expand_decomp.o fortran_file.o fortran_file_vec.o program_line.o utils.o
PROG=kpp4palm.exe

.C.o:
	$(CXX) -c $<

$(PROG): main.o $(OBJ)
	$(CXX) -o $(PROG) main.o $(OBJ)


utils.o: utils.h utils.C
program_line.o:  utils.o program_line.h program_line.C
fortran_file.o: program_line.o fortran_file.h utils.h fortran_file.C
fortran_file_vec.o: program_line.o fortran_file.h utils.h fortran_file_vec.C
expand_decomp.o: fortran_file.o expand_decomp.h expand_decomp.C
create_kpp_module.o: expand_decomp.o create_kpp_module.h create_kpp_module.C
main.o: create_kpp_module.o main.C


clean:
	rm -f *.o $(PROG)
