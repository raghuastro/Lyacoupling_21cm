
!! SED for pop III
!!
subroutine PLANCK_Bnu(nu,T, res)
	use param
	IMPLICIT NONE
	REAL(8),INTENT(IN)::nu,T
	REAL(8)::res
	res=(2.d0*hplanck_ev*nu*nu*nu/C_light/C_light)/(EXP(hplanck_ev*nu/Kboltz_ev/T)-1.d0)

END SUBROUTINE PLANCK_Bnu



!!This will create the SED from the input parameters, for example it will normalize the SED using fx and alpha. 
!!This SED will be used for 1D radiative trasfer. The arrays of SED used for RT should be ev, lum (ev/s/ev)
!!For writtonf the SED for plotting: Must be in the format lamda (A), Intensity (ev/s/lamda)
!
SUBROUTINE generate_SED(n_sed, lambda_s, lum_s)
	use param
	IMPLICIT NONE
	integer, intent(in) :: n_sed
	real(8),dimension(n_sed), intent(out)::lum_s, lambda_s	
	integer::i, j
	real(8),dimension(n_sed)::eng_s, dellambda
	REAL(8)::xz,met,eng, nu, Teff, delnu, tot, norm_source, res
	character(len=15)::s0
	character(180) :: ifile
	real(8), parameter :: sed_factor = 0.95

!! sed for garrelt black body spectrum

		eng=25.d0	!ev
		Teff=3e4
		delnu=1.d0/hplanck_ev	!!delnu for del_E=1ev
		write(*,*) 'Teff=',  Teff
		DO i=1,n_sed
			eng=eng*sed_factor
			eng_s(i)=eng
			nu=eng/hplanck_ev

			call PLANCK_Bnu(nu,Teff, res)

			lum_s(i)=res*delnu
			if(lum_s(i) .lt. 1d-100) then	
				DO j=i+1,n_sed
					eng=eng*sed_factor
					eng_s(j)=eng
				END DO
				lum_s(i+1:n_sed)=0.d0
				exit
			END IF
		END DO


		tot=0.d0
		DO i=1, n_sed-1
		if(eng_s(i)>E_lyalpha) then
		tot=tot+lum_s(i)*(eng_s(i)-eng_s(i+1))
		end if
		END DO

		lum_s=lum_s*5.87d-3*g_gamma_c2ray*L_solar_erg*erg_ev/tot	!!This correspond to a one solar mass halo which emit 1.33d43 number of photon as in sumans paper.
		write(*,*) 'maxval E, L=', maxval(eng_s), maxval(lum_s)

		lambda_s=hplanck_ev*c_light/eng_s*1d8	!!in angstrom
		dellambda=hplanck_ev*c_light*1d8/eng_s/eng_s	!!del lambda for 1eV

		lum_s=lum_s/dellambda/erg_ev	!!erg/A/s
		lum_s=lum_s/f_esc_hi	!! This is the intrinsic SED 


		ifile = output_path//'stellar.dat'
		call write_sed(lambda_s,lum_s,ifile, n_sed)



END SUBROUTINE generate_SED



!!Write the SED of the source
!
SUBROUTINE write_sed(lambda_s,lum_s,filename, n)
	IMPLICIT NONE
	integer, intent(in)::n
	real(8), dimension(n), intent(in)::lambda_s,lum_s
	CHARACTER(len=*) :: filename

	integer::i
	real(8), dimension(n,2)::arrsed
	logical, parameter::plot_sed=.true.


	if(plot_sed) then
		open(unit=41,file=TRIM(filename), status='replace', form='formatted')
		write(*,*) 'SED will be written in', TRIM(filename)
		do i=1,n
			write(41,*) lambda_s(i), lum_s(i)
		end do
		close(41)
	end if

!stop
END SUBROUTINE write_sed




!! This subroutine calculates the Lyman-alpha coupling coefficient.
!! This takes the halos subset (myoffset,myonset) from rank 'myid'.
! 
subroutine visual_calpha(myoffset,myonset,myid)
	use param
	use func
	IMPLICIT NONE
	INTEGER ::myoffset,myonset,myid
	real(8),dimension(nsed)::lamda,lum
	real(8),dimension(n_cell)::arr
	integer,dimension(3)::x1,r
	real(8)::xz,sim_unit,max_dist,res,vol1,r_dist, volmin
	integer::i,j,k,no1,no2,pq,num,q, r1, r2, r3, imin, imax
  	character(180) :: ifile,ofile


	if((myoffset >0) .and. (myonset>=myoffset)) then
!write(*,*) 'myoffset, myonset',myoffset, myonset

		!call zeffectlyman(comm_redshift,arr)
		call zeffectlyman(lya_redshift,arr)
		write(*,*) 'visual_calpha called.. nalpha=',maxval(arr), arr(1), arr(2)

		!sim_unit=box/dble(n_cell)/hlittle/(1.0+comm_redshift) !mpc physical
		sim_unit=box/n_cell/hlittle/(1.0+lya_redshift) !mpc physical

		imax = imax_comm
		imin = imin_comm


		volmin=4.0/3.0*pi*dble(imin)**3.0


		write(*,*) 'imin, imax, time_since, time_since-gap_timestep,',imin, imax, time_since, time_since-gap_timestep

		matrix_alpha=0.0!1d-10

		DO i=myoffset,myonset  !halo
			x1(1:3)=int(dat_overlap(2:4,i)) + 1

		no1=int(volmin)
		if(no1<1) no1=1
		if(imin<1) imin=1

			DO j=imin,imax!int(max_dist)
				r_dist=dble(j)
				res=arr(j)*dble(dat_overlap(1,i))
				if(res < 1d-4) exit
				vol1=4.0/3.0*pi*r_dist**3.0
				no2=int(vol1)
				if(no2>max_cell_read) exit

				IF(no2 >= no1) then
					DO q=no1,no2
						r(:)=x1(:)+cube(q,:)
						call check_boundary(r, n_cell, r1, r2, r3)
						matrix_alpha(r1,r2,r3)=matrix_alpha(r1,r2,r3) + real(res)
					END DO
					no1=no2+1	
				END IF
			END DO
!write(*,*) 'halocalpha calculated sum', i, sum(matrix_alpha), dat_overlap(1,i)
		END DO !halo
	END IF

	write(*,*) 'max,min, sum alpha',maxval(matrix_alpha(:,:,:)),minval(matrix_alpha), sum(matrix_alpha)
!stop

end subroutine visual_calpha

!! This subruitne will call the SED per solar mass and calculate the number of Lyman alpha photon 
!! This uses 1/r/r relation to calculate x_alpha
!!
SUBROUTINE zeffectlyman(z,alpha)
	use param
	use func
	IMPLICIT NONE
	real(8),intent(in)::z
	real(8),dimension(n_cell),intent(out)::alpha
	real(8),dimension(nsed)::lamda,lum
	real(8)::r,zr,lr,l1,l0,e1,e0,en,nu_al,del_l,rp
	integer::i,j,k
	real(8)::xz,res,del_nu, nlya_emission


	lamda = lambdaarr
	lum = lumarr

	del_nu=c_light/angstrom_cm/lamda_alpha/lamda_alpha !!frequency difference for del lamda = 1 A
	nu_al=E_lyalpha/hplanck_ev !Lyman alpha freq in Hz

	alpha=0.0


	DO i=1,n_cell	!!Do we need to increase the size.. photon escaping?? One can try using a increased value
		r=dble(i)*box/hlittle/dble(n_cell) !comoving distace at grid points in Mpc
		zr=rtoz(r,z)
		lr=lamda_alpha*(1.0+zr)/(1.0+z)	!!wavelength in the continuum SEd which redshifted to Lyman alpha at distance r
		if(lr<lamda_hi) then
			alpha(i) = 0.0	!!above EHI, photons will be absorbed for ionization
			write(*,*) 'HI level reached at grid', i
			exit
		else 
			DO j=1,nsed
				if(lamda(j)>=lr) then
					l1=lamda(j)
					l0=lamda(j-1)
					e1=lum(j)
					e0=lum(j-1)
					exit
				end if
			END DO

			en=e0+ (lr-l0)*(e1-e0)/(l1-l0) !linear interpolation !energy per A per sec in erg
			del_l=(1.0+zr)/(1.0+z)! wave length differ for which  at r the wavelength differ will be 1 A!c_light/nu_al/nu_al*(1+z)/(1+zr)*1e8 !del lamda in A for unit frequency at r
			rp=r/(1.0+zr)*megaparsec
			alpha(i)=en*del_l/4.0/pi/rp/rp/del_nu*erg_ev/E_lyalpha
			alpha(i)=1.66e11/(1.0+zr)*alpha(i)*f_esc_lya*(1.0+0.6*(1.0-f_esc_hi))  ! ly coeffi for unit mass

		end if
	END DO




END SUBROUTINE zeffectlyman

