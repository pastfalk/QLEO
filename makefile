objects=param_mod.o main.o get_deltaf.o splcoeff_spec.o Term1.o Term2.o Term3.o Term4.o Term5.o Term6.o Term1_pole.o Term2_pole.o Term3_pole.o Term4_pole.o Term5_pole.o Term6_pole.o cubic_sol.o read_dist.o deriv_dist.o read_data.o disp_rel.o polyfit.o muller.o disp_det.o exp_Bessel_In.o exp_Bessel_In_mpfun.o exp_dBessel_In.o gamma_func.o Z_func.o dZ_func.o integrator.o splcoeff_dist.o spline_interpol.o Bessel_int.o acc_F.o F23.o F12.o F23_mpfun.o F12_mpfun.o acc_Kvpa.o int_para.o int_para_mpfun.o cerror.o cont_frac.o splcoeff_beta.o get_beta.o fort_Bes.o print_dist.o print_omega.o print_df.o initiate.o compute_EB.o adapt_krange.o
mpfun_obj = mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq2.o  mpmodule.o
mpfun_mod = mpfuna.mod  mpfunb.mod mpfunc.mod mpfund.mod mpfune.mod mpfunf.mod mpfung.mod  mpmodule.mod
f90comp = gfortran
options =  -fdefault-real-8 -g -O3 -fopenmp -fbounds-check -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero #-ffpe-trap=denormal -Wall
options_mp = -O3 -g #-ffpe-trap=invalid -ffpe-trap=zero -ffpe-trap=overflow #-ffpe-trap=denormal

qsolve: $(objects) $(mpfun_obj)
	$(f90comp) -o qsolve $(options) $(objects) $(mpfun_obj)

param_mod.mod: param_mod.o param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

param_mod.o: param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

mpmodule.mod: mpmodule.o mpmodule.f90
	$(f90comp) -c $(options_mp) mpmodule.f90

mpmodule.o: mpmodule.f90
	$(f90comp) -c $(options_mp) mpmodule.f90

mpfuna.mod: mpfuna.o mpfuna.f90
	$(f90comp) -c $(options_mp) mpfuna.f90

mpfuna.o: mpfuna.f90
	$(f90comp) -c $(options_mp) mpfuna.f90

mpfunb.mod: mpfunbq.o mpfunbq.f90
	$(f90comp) -c $(options_mp) mpfunbq.f90

mpfunbq.o: mpfunbq.f90
	$(f90comp) -c $(options_mp) mpfunbq.f90

mpfunc.mod: mpfunc.o mpfunc.f90
	$(f90comp) -c $(options_mp) mpfunc.f90

mpfunc.o: mpfunc.f90
	$(f90comp) -c $(options_mp) mpfunc.f90

mpfund.mod: mpfund.o mpfund.f90
	$(f90comp) -c $(options_mp) mpfund.f90

mpfund.o: mpfund.f90
	$(f90comp) -c $(options_mp) mpfund.f90

mpfune.mod: mpfune.o mpfune.f90
	$(f90comp) -c $(options_mp) mpfune.f90

mpfune.o: mpfune.f90
	$(f90comp) -c $(options_mp) mpfune.f90

mpfunf.mod: mpfunf.o mpfunf.f90
	$(f90comp) -c $(options_mp) mpfunf.f90

mpfunf.o: mpfunf.f90
	$(f90comp) -c $(options_mp) mpfunf.f90

mpfung.mod: mpfungq2.o mpfungq2.f90
	$(f90comp) -c $(options_mp) mpfungq2.f90

mpfungq2.o: mpfungq2.f90
	$(f90comp) -c $(options_mp) mpfungq2.f90

main.o: param_mod.mod main.f90
	$(f90comp) -c $(options) main.f90

get_deltaf.o: param_mod.mod get_deltaf.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal get_deltaf.f90

splcoeff_spec.o: param_mod.mod splcoeff_spec.f90
	$(f90comp) -c $(options) splcoeff_spec.f90

Term1.o: param_mod.mod Term1.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term1.f90

Term2.o: param_mod.mod Term2.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term2.f90

Term3.o: param_mod.mod Term3.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term3.f90

Term4.o: param_mod.mod Term4.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term4.f90

Term5.o: param_mod.mod Term5.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term5.f90

Term6.o: param_mod.mod Term6.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term6.f90

Term1_pole.o: param_mod.mod Term1_pole.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term1_pole.f90

Term2_pole.o: param_mod.mod Term2_pole.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term2_pole.f90

Term3_pole.o: param_mod.mod Term3_pole.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term3_pole.f90

Term4_pole.o: param_mod.mod Term4_pole.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term4_pole.f90

Term5_pole.o: param_mod.mod Term5_pole.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term5_pole.f90

Term6_pole.o: param_mod.mod Term6_pole.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal Term6_pole.f90

cubic_sol.o: param_mod.mod cubic_sol.f90
	$(f90comp) -c $(options) cubic_sol.f90

read_dist.o: param_mod.mod read_dist.f90
	$(f90comp) -c $(options) read_dist.f90

read_data.o: param_mod.mod read_data.f90
	$(f90comp) -c $(options) read_data.f90

deriv_dist.o: param_mod.mod deriv_dist.f90
	$(f90comp) -c $(options) -ffpe-trap=invalid -ffpe-trap=overflow -ffpe-trap=zero -ffpe-trap=denormal deriv_dist.f90

disp_rel.o: param_mod.mod disp_rel.f90
	$(f90comp) -c $(options) disp_rel.f90

polyfit.o: polyfit.f90
	$(f90comp) -c $(options) polyfit.f90

muller.o: param_mod.mod muller.f90
	$(f90comp) -c $(options) muller.f90

disp_det.o: param_mod.mod disp_det.f90
	$(f90comp) -c $(options) disp_det.f90

exp_Bessel_In.o: param_mod.mod exp_Bessel_In.f90
	$(f90comp) -c $(options) exp_Bessel_In.f90

exp_Bessel_In_mpfun.o: param_mod.mod  $(mpfun_mod) exp_Bessel_In_mpfun.f90
	$(f90comp) -c $(options)  exp_Bessel_In_mpfun.f90

exp_dBessel_In.o: exp_dBessel_In.f90
	$(f90comp) -c $(options) exp_dBessel_In.f90

gamma_func.o: param_mod.mod gamma_func.f90
	$(f90comp) -c $(options) gamma_func.f90

Z_func.o: param_mod.mod Z_func.f90
	$(f90comp) -c $(options) Z_func.f90

dZ_func.o: dZ_func.f90
	$(f90comp) -c $(options) dZ_func.f90

integrator.o: param_mod.mod integrator.f90
	$(f90comp) -c $(options) integrator.f90

splcoeff_dist.o: param_mod.mod splcoeff_dist.f90
	$(f90comp) -c $(options) splcoeff_dist.f90

spline_interpol.o: spline_interpol.f90
	$(f90comp) -c $(options) spline_interpol.f90

Bessel_int.o: param_mod.mod Bessel_int.f90
	$(f90comp) -c $(options) Bessel_int.f90

acc_F.o: acc_F.f90
	$(f90comp) -c $(options) acc_F.f90

acc_Kvpa.o: acc_Kvpa.f90
	$(f90comp) -c $(options) acc_Kvpa.f90

F23.o: F23.f90
	$(f90comp) -c $(options) F23.f90

F12.o: F12.f90
	$(f90comp) -c $(options) F12.f90

F23_mpfun.o: $(mpfun_mod) F23_mpfun.f90
	$(f90comp) -c $(options) F23_mpfun.f90

F12_mpfun.o: $(mpfun_mod) F12_mpfun.f90
	$(f90comp) -c $(options) F12_mpfun.f90

int_para.o: param_mod.mod int_para.f90
	$(f90comp) -c $(options) int_para.f90

int_para_mpfun.o: param_mod.mod  $(mpfun_mod) int_para_mpfun.f90
	$(f90comp) -c $(options) int_para_mpfun.f90

cerror.o: param_mod.mod cerror.f90
	$(f90comp) -c $(options) cerror.f90

cont_frac.o: cont_frac.f90
	$(f90comp) -c $(options)  cont_frac.f90

get_beta.o: param_mod.mod get_beta.f90
	$(f90comp) -c $(options) get_beta.f90

splcoeff_beta.o: param_mod.mod splcoeff_beta.f90
	$(f90comp) -c $(options) splcoeff_beta.f90

fort_Bes.o: fort_Bes.f90
	$(f90comp) -c $(options)  fort_Bes.f90

print_dist.o: param_mod.mod print_dist.f90
	$(f90comp) -c $(options) print_dist.f90

print_df.o: param_mod.mod print_df.f90
	$(f90comp) -c $(options) print_df.f90

print_omega.o: param_mod.mod print_omega.f90
	$(f90comp) -c $(options) print_omega.f90

initiate.o: param_mod.mod initiate.f90
	$(f90comp) -c $(options) initiate.f90

compute_EB.o: param_mod.mod compute_EB.f90
	$(f90comp) -c $(options) compute_EB.f90

adapt_krange.o: param_mod.mod adapt_krange.f90
	$(f90comp) -c $(options) adapt_krange.f90

clean:
	rm $(mpfun_obj)
	rm $(mpfun_mod)
	rm param_mod.mod
	rm $(objects)
	rm qsolve 
