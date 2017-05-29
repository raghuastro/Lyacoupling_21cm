subroutine count_num_halo(myid, z, nhalo)
	use param
	use param_halo
	implicit none
	integer::myid
	real(8)::z

  	character(7) :: z_s
  	character(250) :: ifile
	real(8) ::halo_max,halo_min,M,Mass, Meff, mass_gal, M1, z_ns
	integer::nhalo,i,j,k,cnt,nhalo_max, ll, xx, yy, cnt1, fstat, nhalo_ns

	integer, dimension(n_cell, n_cell, n_cell)::halo_trace_ns
	real(8), dimension(n_cell, n_cell, n_cell)::halo_mass_new_ns

	halo_trace=0 !! This tracks the grid with halo..1 mean halo present

!!choose the halo file : consider feedback, large/small source

if(test_signle_source) then
nhalo=1
else

	
	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=halolist_path//z_s(1:len_trim(z_s))//"-coarsest_sources.dat"
	write(*,*) 'halo file, z', ifile, z
	open(unit=13,file=ifile, status='old')
	read(13, *) nhalo
	close(13)
end if

end subroutine count_num_halo


!!This subroutine reads the list of halo at redshift z

!!This track age
!!The halo file should be given in grid positions, sum all haloes >Mcut at every grid points and place.
!!Sould not have many similar sources at the same grid points, otherwise age calculation will be wrong.
!!
subroutine read_halofile_track_age(myid, z)
	use param
	use param_halo
	implicit none
	integer::myid
	real(8)::z

  	character(7) :: z_s
  	character(250) :: ifile
	real(8) ::halo_max,halo_min,M,Mass, Meff, mass_gal, M1, z_ns
	integer::nhalo,i,j,k,cnt,nhalo_max, ll, xx, yy, cnt1, fstat, nhalo_ns

	integer, dimension(n_cell, n_cell, n_cell)::halo_trace_ns
	real(8), dimension(n_cell, n_cell, n_cell)::halo_mass_new_ns

	real(8) :: m0, gama, zi

	halo_trace=0 !! This tracks the grid with halo..1 mean halo present

if(test_signle_source) then
	gama=0.0
	!gama=0.5
	m0=1d11
	zi=25.0

	mass=m0*exp(gama*(zi-z))


!	Mhalo_tot=Mhalo_tot+mass

	num_halo=1
	cnt=1
	i=n_cell/2
	j=i
	k=i

	halo_mass_new=0.0
	halo_trace(i,j,k)=1

	Meff = mass
	halo_mass_new(i,j,k)=Meff
	write(*,*) 'effective mass, normal mass, age', Meff, mass, halo_age_new(i,j,k)
	dat_overlap(1,cnt)=Meff !!mass
	dat_overlap(2,cnt)=real(i)+0.5	
	dat_overlap(3,cnt)=real(j)+0.5	
	dat_overlap(4,cnt)=real(k)+0.5



else

!!choose the halo file : consider feedback, large/small source

	
	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=halolist_path//z_s(1:len_trim(z_s))//"-coarsest_sources.dat"
	write(*,*) 'halo file', ifile
	open(unit=13,file=ifile, status='old')
	read(13, *) nhalo_max

	write(*,*) 'Max number of halo', nhalo_max


	nhalo=0
	halo_max=0.0
	halo_min=1.0d+20
	halo_mass_new=0.0		!!This is the new array of halo mass


!!Read the haloes :

	DO ll=1, nhalo_max ! halo loop
		read(13,*)  i,j,k,M,M1
		mass=M*M_grid_s  ! mass of the halo in solar mass

	if(mass >1d8) then
		halo_trace(i,j,k)=1 !To identify the cells with source	
		mass_gal=Mass
			
		nhalo=nhalo+1
		if(mass_gal .lt. halo_min) halo_min=mass_gal
		if(mass_gal .gt. halo_max) halo_max=mass_gal


		halo_mass_new(i,j,k)=halo_mass_new(i,j,k)+real(mass_gal)

	end if
				
	END DO !l
	print*,'Halo considered',nhalo


	close(13)
	num_halo=nhalo

	write(*,*) 'mim max mass halo', halo_min, halo_max
!stop



!!Now we will calculate the effective mass of the halo over time, this is to account the mass evlution of the halo which is in general exponential.
!!However, we will not taken into account the dennisty contrast evolution around the source over the time and will work only with the current density distribution at that redshift.

	dat_overlap=0.0	!!This array will carry the mass, X(:), etc.. 
	cnt=0


	DO k=1,n_cell
	DO j=1,n_Cell
	DO i=1,n_cell
		xx=halo_trace(i,j,k)
		if(xx ==1 ) then
			cnt=cnt+1
			Meff=halo_mass_new(i,j,k)

			dat_overlap(1,cnt)=Meff 
			dat_overlap(2,cnt)=real(i)+0.5	
			dat_overlap(3,cnt)=real(j)+0.5	
			dat_overlap(4,cnt)=real(k)+0.5
		end if
	END DO
	END DO
	END DO

	print*,'Read_halofile exit. nhalo',cnt

end if


end subroutine read_halofile_track_age


