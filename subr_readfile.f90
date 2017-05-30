!!!!read the halofinds redshifts for the simulation
!
subroutine read_redshift
use param
	implicit none
	integer::fstat,i
	real(4)::arg


	open(11,file=checkpoints,status='old',iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!ERROR!!', checkpoints
		stop
	end if

	DO num_checkpoints=1,max_input
	    read(unit=11,err=51,end=41,fmt='(f20.10)') arg
		z_checkpoint(num_checkpoints)=dble(arg)
	END DO
	41  num_checkpoints=num_checkpoints-1
	51  close(11)


	write(*,*) 'halofinds to recompose:'
	do i=1,num_checkpoints
	    write(*,'(f5.1)') z_checkpoint(i)
	enddo



end subroutine read_redshift

!!calculate the time between the simulation snapshots..
!!
subroutine cal_age_timestep
	use param
	implicit none
	real(8)::z1, z2, delt
	integer::i

	z2=z_checkpoint(1)

	DO i=1,num_checkpoints-1
		z1=z_checkpoint(i+1)
		call timegap(z1,z2,delt)
		age_checkpoint(i)=delt
	END DO
	age_checkpoint(num_checkpoints)=age_checkpoint(num_checkpoints-1) + 10.0

	write(*,*) 'Age to recompose:'
	do i=1,num_checkpoints
	    write(*,'(f5.1)') age_checkpoint(i)
	enddo

end subroutine cal_age_timestep


!!This will generate an array of nearby grid points from 0,0,0 within a box of dimension n,n,n
!!
SUBROUTINE readcube_new(n)
	use param
	implicit none
	integer, intent(in) ::n
	integer :: nin
	real(4), dimension(:), allocatable::arr, arr1
	integer, dimension(:,:), allocatable ::cube1
	integer::i,j,k,ii,jj,kk,cnt, stat, cn, fstat
	real(4) ::d, tmp
	integer, dimension(3)::temp
	character(180) :: ifile
	character(12) :: grid_s

	nin=n/2
	allocate(arr(nin*nin*nin), arr1(nin*nin*nin), cube1(nin*nin*nin,3))

	ii=1
	jj=1
	kk=1
	cnt=0

	DO i=1, nin
	DO j=1, nin
	DO k=1, nin
		cnt=cnt+1
		d=real(i-ii)**2.0 + real(j-jj)**2.0 + real(k-kk)**2.0
		d=sqrt(d)
		arr(cnt) = d

		cube1(cnt,1) = i-ii
		cube1(cnt,2) = j-jj
		cube1(cnt,3) = k-kk
	END DO
	END DO
	END DO

	arr1=arr


	write(*,*) 'Before'
	write(*,*) arr(1), arr(2), maxval(arr)

	call sort(nin*nin*nin,arr)

	write(*,*) 'After'
	write(*,*) arr(1), arr(2), arr(cnt)

	cn=0
	DO i=1, cnt
		!if(mod(i,10000)==0) write(*,*) i
		d=arr(i)
		temp=cube1(i,:)
		tmp=arr1(i)
		stat=0
		DO j=i,cnt
		if(d==arr1(j)) then
			ii=cube1(j,1)
			jj=cube1(j,2)
			kk=cube1(j,3)
			stat=1
			arr1(j) =1e30
			arr1(i)=arr1(j)
			arr1(j)=tmp
			cube1(j,:)=temp


			cn=cn+1
			cube(cn,:)=(/-ii,-jj,-kk/)

			if(kk .ne. -kk) then
				cn=cn+1
				cube(cn,:)=(/-ii,-jj,kk/)
			end if
			if(jj .ne. -jj) then
				cn=cn+1
				cube(cn,:)=(/-ii,jj,-kk/)
			end if
			if(ii .ne. -ii) then
				cn=cn+1
				cube(cn,:)=(/ii,-jj,-kk/)
			end if
			if(jj .ne.  -jj .and. kk .ne. -kk) then
				cn=cn+1
				cube(cn,:)=(/-ii,jj,kk/)
			end if

			if(ii .ne.  -ii .and. kk .ne. -kk) then
				cn=cn+1
				cube(cn,:)=(/ii,-jj,kk/)
			end if
			if(jj .ne.  -jj .and. ii .ne. -ii) then
				cn=cn+1
				cube(cn,:)=(/ii,jj,-kk/)
			end if

			if(jj .ne.  -jj .and. kk .ne. -kk .and. ii .ne. -ii) then
				cn=cn+1
				cube(cn,:)=(/ii,jj,kk/)
			end if

			exit
		end if
		END DO
		if(stat==0) then
		write(*,*) 'not convergeing'
		stop
		end if
	END DO

	write(grid_s,'(I3)') ncube
	grid_s=adjustl(grid_s)


	ifile=cellpath//'cube'//trim(grid_s)//'.bi'
	write(*,*) 'cube array writting to file', ifile
	open(unit=311,file=ifile, form='binary', status='replace', iostat=fstat)
	write(311) cn
	write(311) cube
	close(311) 

	write(*,*) 'cube array dim', cn, n*n*n, n, real(cn)**(1./3.)

	write(*,*) cube(1,:), cube(2,:), cube(cn,:)
	deallocate(arr, arr1, cube1)

END SUBROUTINE readcube_new


!! This reads the array of nearby points from 0,0,0.. in case file misiing will generate the array
!!
SUBROUTINE readcube
	use param
	implicit none
  	character(180) ::ofile5,ifile1, ifile
	integer::q,q1,count,i,j,k,fs,n2p,count1,num, fstat
	character(12) :: grid_s

	write(grid_s,'(I3)') ncube
	grid_s=adjustl(grid_s)


	ifile=cellpath//'cube'//trim(grid_s)//'.bi'
	open(unit=31,file=ifile, form='binary', status='old', iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!ERROR!', ifile
		write(*,*) 'File will be regenerated'
		call	readcube_new(ncube)
	else
		read(31) q
		if(q .ne. max_cell_read) then
		write(*,*) 'Dimension of cube array not matching.. terminating job.', q, max_cell_read
		stop
		end if
		read(31) cube
	end if	

					
	close(31)
	write(*,*) 'done cube reading, maxval cube',maxval(cube), cube(1,:),  cube(2,:)


END SUBROUTINE readcube

!!To check periodicity in the box
!!
subroutine check_boundary(x, n, i, j, k)
	implicit none

	integer, dimension(3), intent(in)::x
	integer, intent(in)::n
	integer, intent(out)::i,j,k

	i=x(1)
	j=x(2)
	k=x(3)

	if(i<1) i=i+n
	if(j<1) j=j+n
	if(k<1) k=k+n

	if(i>n) i=i-n
	if(j>n) j=j-n
	if(k>n) k=k-n

end subroutine check_boundary


