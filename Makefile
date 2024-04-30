OUTPUT_NAME=QMC_run
FC = gfortran

$(OUTPUT_NAME): QMC_main.F90 QMC_energy.F90 QMC_utilities.F90
	$(FC) -o $@ $^

clean:
	rm -f $(OUTPUT_NAME)

