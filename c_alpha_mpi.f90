!! make -j8 -f Make_c_alpha
!! mpirun -np 4 ./c_alpha

!!This program calculate the x_alpha maps for a given source distribution
!!Assumption 1/r^2 fall of the Lyalpha photns number density
!!
program mpi_sim
	use param
	use param_halo
	implicit none 
     	include 'mpif.h'
      	integer  ierr, offset, onset, i, j, k, tag1, &
     	&         tag2, tag3, tag4, source, check,count,check_min
      	real(8) ::   mysum, mysum1, total_sum,total_sum1,af, sim_unit, m_tot
	real(4) :: time
      	integer  status(MPI_STATUS_SIZE),fstat
  	character(180) :: ifile,ofile
  	character(7) :: z_s, id_s
	logical,parameter::restart=.false.		!! To restart at some particular slice
	INTEGER,PARAMETER::restart_id=15
	real(4), dimension(n_cell, n_cell, n_cell):: dummy_arr


!C ***** Initializations *****
      	call MPI_INIT(ierr)
      	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
      	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

		tag1=1
		tag2=2
		tag3=3
		tag4=4


!C***** Master task only ****** read
      if (rank .eq. 0) then
		call readcube
		call read_redshift
		call cal_age_timestep

		call generate_SED(nsed, lambdaarr, lumarr)
		halo_age_new=0
		halo_mass_old=0.0
		write(*,*) '# of cell =', n_cell
		write(*,*) '# of checkpoints = ',num_checkpoints
	
		open(unit=16, file=output_path//'x_alpha.dat')


		call MPI_BCAST(num_checkpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(cube,max_cell_read*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(lambdaarr,nsed,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(lumarr,nsed,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	else
		call MPI_BCAST(num_checkpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(cube,max_cell_read*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(lambdaarr,nsed,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(lumarr,nsed,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	end if

	if(restart) then
		write(*,*) 'Program is restarted : redshift',z_checkpoint(restart_id)
		check_min=restart_id
	else
		check_min=1
	end if

DO check=1, num_checkpoints

	if(rank==0) then
		write(*,*) '****************'
		write(*,*) 'Processing redshift = ',z_checkpoint(check)
		lya_redshift=z_checkpoint(check)	!!This is the redshift where the coupling will be calculated..
		matrix_alpha=0.0!1d-10
		dummy_arr=0.0!1d-10
		sim_unit=box/n_cell/hlittle/(1.0+lya_redshift) !mpc physical
	end if

!!start the loop to consider the cntributions from previous redshifts..

DO redshift_loop_index=1,check

	if(rank==0) then
		comm_redshift=z_checkpoint(redshift_loop_index)
		write(*,*) 'z loop start, index, z', redshift_loop_index, comm_redshift
		gap_timestep=age_checkpoint(1)
		time_since= age_checkpoint(check)-age_checkpoint(redshift_loop_index)+gap_timestep	!time gap between zstar and zstep


		if(redshift_loop_index==1)	then
			imax_comm = int((time_since)*Megayear*c_light/Megaparsec/sim_unit)!*(1.0+lya_redshift)/(1.0+comm_redshift)
		else
			imax_comm = imin_comm
		end if

		imin_comm = int((time_since-gap_timestep)*Megayear*c_light/Megaparsec/sim_unit)!*(1.0+lya_redshift)/(1.0+comm_redshift)
		write(*,*) 'timesince, timegap', time_since, gap_timestep
		write(*,*) 'imin, imax, i, check', imin_comm, imax_comm, i, check

		call count_num_halo(rank, comm_redshift, num_halo)

		call MPI_BCAST(num_halo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(comm_redshift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(lya_redshift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(time_since,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(gap_timestep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(imax_comm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(imin_comm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	else
		call MPI_BCAST(num_halo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(comm_redshift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(lya_redshift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(time_since,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(gap_timestep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(imax_comm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(imin_comm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	end if


	IF(rank ==0) THEN
	if(num_halo>0) then
		allocate(dat_overlap(4,num_halo))
		call read_halofile_track_age(rank, comm_redshift)
		call cpu_time(time)
		write(*,*) 'Halofile reading : done : time, num_halo ',time, num_halo


		DO i=1,numtasks-1
			call para_range(num_halo, numtasks, i, offset, onset)
		     	call MPI_SEND(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
		     	call MPI_SEND(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr)
			if((offset >0) .and. (onset>=offset)) then
			     	call MPI_SEND(dat_overlap(1,offset),4*(onset-offset+1), MPI_REAL, i, tag4, MPI_COMM_WORLD, ierr) 
			     	call MPI_SEND(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr) 
			end if
		end do 


		call para_range(num_halo, numtasks, rank, offset, onset)
		call visual_calpha(offset,onset,rank)


		write(*,*) 'before:sum alpha this loop', sum(dummy_arr)

!if(redshift_loop_index==1) then

		dummy_arr=dummy_arr+matrix_alpha
		call cpu_time(time)
		write(*,*) 'RANK0 : done :time ', time
!end if

		do i=1, numtasks-1
		     	call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
		     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
			if((offset >0) .and. (onset>=offset)) then
		     	call MPI_RECV(dat_overlap(1,offset),4*(onset-offset+1), MPI_REAL, i, tag4, MPI_COMM_WORLD, status, ierr)
		     	call MPI_RECV(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
			dummy_arr=dummy_arr+matrix_alpha
			end if
			write(*,*) 'Recv,rank=',i
		end do 

		m_tot = sum(dat_overlap(1,1:num_halo))

		write(*,*) 'after:sum alpha this loop', sum(dummy_arr)
		deallocate(dat_overlap)

		write(id_s,'(f7.3)') comm_redshift
		id_s=adjustl(id_s)

		write(z_s,'(f7.3)') z_checkpoint(check)
		z_s=adjustl(z_s)


if(test_signle_source) then
!		ofile=output_path//trim(z_s)//"matrix_alphastar_slice.dat"//trim(id_s)
!		write(*,*) ofile
!		call write_slice(dummy_arr, ofile, n_cell, n_cell, n_cell, n_cell/2)
end if

	if(redshift_loop_index==check) then
	
		matrix_alpha=dummy_arr
		write(*,*) 'max,min alpha',maxval(matrix_alpha(:,:,:)),minval(matrix_alpha)
		count=0
		DO i=1,n_cell
		DO j=1,n_cell
		DO k=1,n_cell
			if(matrix_alpha(i,j,k)<10.0) then
			count=count+1
			end if
		END DO
		END DO
		END DO

		af=1.0-dble(count)/dble(n_cell)**3.0
		write(*,*) 'coupled frac=',af


		ofile=output_path//trim(z_s)//"matrix_alphastar.dat"
		write(*,*) ofile
		call write_binary_4bit(matrix_alpha, ofile, n_cell, n_cell, n_cell)

		if(test_signle_source) then
			ofile=output_path//trim(z_s)//'1dalpha.dat'
			write(*,*) ofile
			open(unit=24, file=ofile, form='formatted')

			DO i=n_cell/2, n_cell
			write(24,'(17e12.4)') dble(i-n_cell/2), matrix_alpha(i,n_cell/2, n_cell/2)
			END DO
			 close(24)
		end if
!
!		ofile=output_path//trim(z_s)//"matrix_alphastar_slice.dat"
!		write(*,*) ofile
!		call write_slice(matrix_alpha, ofile, n_cell, n_cell, n_cell, n_cell/2)

		call cpu_time(time)
		write(*,*) 'Finish : time, z', time, comm_redshift

		write(16,'(17e12.4)') comm_redshift, af, m_tot, sum(matrix_alpha)

		if(af> 0.99) then
			write(*,*) 'coupling completed'
			close(16)
			call mpi_abort(mpi_comm_world, ierr)
		end if
	end if	!!
	end if
	ELSE
		if(num_halo>0) then
		allocate(dat_overlap(4,num_halo))
		end if

	END IF 

 
 

	if(rank .ne. 0) then
		call MPI_RECV(offset, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(onset, 1, MPI_INTEGER, 0, tag2, MPI_COMM_WORLD, status, ierr)
		if((offset >0) .and. (onset>=offset)) then
		call MPI_RECV(dat_overlap(1,offset),4*(onset-offset+1), MPI_REAL, 0, tag4, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, 0, tag3, MPI_COMM_WORLD, status, ierr)


		call visual_calpha(offset,onset,rank)
		end if


		call MPI_SEND(offset, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, ierr) 
		call MPI_SEND(onset, 1, MPI_INTEGER, 0, tag2, MPI_COMM_WORLD, ierr)
		if((offset >0) .and. (onset>=offset)) then
		call MPI_SEND(dat_overlap(1,offset),4*(onset-offset+1), MPI_REAL, 0, tag4, MPI_COMM_WORLD, ierr)
		call MPI_SEND(matrix_alpha, n_cell*n_cell*n_cell, MPI_REAL, 0, tag3, MPI_COMM_WORLD, ierr)

		end if
		if(num_halo>0) then
		deallocate(dat_overlap)
		end if
	end if

end do

end do

      call MPI_FINALIZE(ierr)


end program


