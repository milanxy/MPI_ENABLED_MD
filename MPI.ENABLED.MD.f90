!**********************************************************************************************!
! A MOLECULAR DYNAMICS PROGRAM ON LJ SYSTEM                                                    !
!**********************************************************************************************!
!**********************************************************************************************!
!MODULE1: PRECISION CONTROL OF THE PROGRAM IS INCORPORATED HERE                                !
!**********************************************************************************************!
      module Precision_Control
      implicit none
      integer, parameter :: single = selected_real_kind(p=6,r=37)    !single precision
      integer, parameter :: double = selected_real_kind(p=15,r=37)  !double precision
      integer, parameter :: quad   = selected_real_kind(p=33,r=4931) !quadruple precision
      save
      end module
!***********************************************************************************************!
!MODULE2: PARAMETERS OF THE PROGRAM                                                             !
!***********************************************************************************************!
      module Parameter_Control
      Use Precision_Control
      Use mpi
      integer::Np ,ndim , Nstep , Nsolvent
      real(kind=double):: density , boxlen , epsil  , init_temp ,temp,dt
      real(kind=double)::rcut
!! variables to use mpi parallel in different nodes
      integer ::ierr , my_id , num_procs , nop
      integer ::root_process , particle_per_process , particle_for_root
      integer ::send_data_tag=2001 , return_data_tag=2002
      integer:: status(MPI_STATUS_SIZE)
      end module
!***********************************************************************************************!
!MODULE3: ALL THE READ    VARIABLES                                                             !
!***********************************************************************************************!
      module Read_Variables
      Use Precision_Control
     
      real(kind=double),allocatable,dimension(:,:):: R ,init_vel , Rm
      end module
!***********************************************************************************************!
!MODULE4: ALL THE   VARIABLES                                                                   !
!***********************************************************************************************!
      module Other_Variables
      Use Precision_Control
      
      real(kind=double),allocatable,dimension(:):: v_sum
      real(kind=double)::v2_sum ,sf
      real(kind =double)::potential
      real(kind=double),allocatable,dimension(:,:)::net_force
      end module
!***********************************************************************************************!
!PROGRAM:                                                                                       !
!***********************************************************************************************!
       program main
       Use Precision_Control
       Use Parameter_Control
       Use Read_Variables
       Use Other_Variables
       use mpi
       implicit none
       integer ::istep , ii 
       !real(kind=double)::
        call INITIALIZE
        call READ_TRAJECTORY
        call VELOCITY_ALLOCATION
       call MPI_INIT ( ierr )
        
       call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
       call MPI_COMM_SIZE (MPI_COMM_WORLD, nop, ierr)
       call FORCE
       !init_vel = 0
      
        
      do istep = 1 , Nstep
         if(my_id .eq. 0) then
         v_sum = 0
         v2_sum = 0
        do ii = 1 , Np
      !     write(200 , *) net_force(ii, 1) , net_force(ii, 2) , net_force(ii, 3)
           R(ii , :) = R(ii , :) + init_vel(ii , :)*dt + 0.50*net_force(ii , :)*dt**2
           if( istep .gt. 320000 ) then
             write(13 ,*)  R(ii , 1) ,  R(ii , 2) ,  R(ii , 3)
           end if

        ! if(istep .gt.1000.and. istep .lt. 20000 ) then
        !   sf = sqrt(init_temp/temp)
        !   init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
        ! else
        !   init_vel(ii , :) = init_vel(ii , :) +0.50*net_force(ii ,:)*dt
        ! end if
        
        if( istep .gt. 1000 .and. istep .le. 20000 ) then
           sf = sqrt(5.0/temp)
           init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
         else if(istep .gt.20000 .and. istep .le. 300000) then
           init_vel(ii , :) = init_vel(ii , :) +0.50*net_force(ii ,:)*dt
         else if(istep .gt. 300000 .and. istep .le. 305000) then
           sf = sqrt(4.0/temp)
           init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
         else if(istep .gt.305000 .and. istep .le. 310000) then
           sf = sqrt(3.0/temp)
           init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
         else if(istep .gt.310000 .and. istep .le. 315000) then
           sf = sqrt(2.0/temp)
           init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
         else if(istep .gt.315000 .and. istep .le. 320000) then
           sf = sqrt(1.0/temp)
           init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
         else if(istep .gt.320000 .and. istep .le. 350000) then
           sf = sqrt(init_temp/temp)
           init_vel(ii , :) = init_vel(ii , :)*sf +0.50*net_force(ii ,:)*dt
         else
           init_vel(ii , :) = init_vel(ii , :) +0.50*net_force(ii ,:)*dt
         end if





        end do
!        end if
        net_force = 0
        end if
        call FORCE
        if(my_id .eq. 0) then
        !  write(40 , *) potential
        do ii = 1 , Np
           init_vel(ii , :) = init_vel(ii , :) +0.50*net_force(ii ,:)*dt
           v_sum(:) = v_sum(:)+init_vel(ii , :)
           !v2_sum = v2_sum +init_vel(ii , 1)**2+ init_vel(ii , 2)**2+init_vel(ii , 3)**2
        end do
          v_sum(:) = v_sum(:)/Real(Np)
        do ii = 1 , Np
           init_vel(ii , :) = (init_vel(ii , :) - v_sum(:))
          ! write(50 , *)  init_vel(ii , 1) , init_vel(ii , 2) , init_vel(ii , 3)
           v2_sum = v2_sum +init_vel(ii , 1)**2+ init_vel(ii , 2)**2+init_vel(ii , 3)**2
        end do
        temp =v2_sum/(3.0*Np)
        !write(11 , *) istep , pot_energy
        write(12 , *) istep , temp
        end if
      end do
        call MPI_FINALIZE(ierr)
 
       end program main
!***********************************************************************************************!
!SUBROUTINE1:INITIALIZATION OF VARIABLES                                                        !
!***********************************************************************************************!
      subroutine INITIALIZE
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Other_Variables
      implicit none
      density   = 0.85
      dt        = 0.002
      init_temp = 0.5
      Np        = 5324
      ndim      = 3
      Nsolvent  =  3194
      boxlen    =  (Np/density)**(0.333333333)
      print*, boxlen
      Nstep     = 2000000
      rcut      = (boxlen/2)**2.0
      
      allocate(R(Np ,ndim) , init_vel(Np , ndim) , Rm(Np , ndim))
      allocate(net_force(Np , ndim))
      allocate(v_sum(3))
      
      R        = 0
      init_vel = 0
      v_sum    = 0
      v2_sum   = 0
      net_force= 0
      potential =0
      temp      = 0
open(unit=1,file="/home/milan/SPINODAL.RUN.CHECK/INITIAL.FCC.LATICE/initial.dat",action="read")
      end subroutine INITIALIZE
!***********************************************************************************************!
!SUBROUTINE2:READING INITIAL COORDINATES                                                        !
!***********************************************************************************************!
      subroutine READ_TRAJECTORY
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Other_Variables
      implicit none
      integer :: ii
      do ii = 1 , Np
         read(1 ,*) R(ii , 1) , R(ii , 2) , R(ii , 3)
      end do
      end subroutine READ_TRAJECTORY
!***********************************************************************************************!
!SUBROUTINE3:MAXWELLIAN VELOCITY ALLOTMENT                                                      !
!***********************************************************************************************!
      subroutine VELOCITY_ALLOCATION
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Other_Variables
      implicit none
      real(kind=double)::a,b,c
      integer::ii
       do ii = 1 ,Np
         call RAN_GAUSS(a,b,c)
         init_vel(ii , 1) = a
         init_vel(ii , 2) = b
         init_vel(ii , 3) = c

       end do
       do ii = 1 ,Np
         v_sum(1) = v_sum(1)+init_vel(ii ,1)
         v_sum(2) = v_sum(2)+init_vel(ii ,2)
         v_sum(3) = v_sum(3)+init_vel(ii ,3)
         v2_sum = v2_sum + init_vel(ii ,1)**2+init_vel(ii ,2)**2+init_vel(ii ,3)**2
       end do
       temp= v2_sum/(3*Np)
       !print*, temp
         v_sum(1) = v_sum(1)/Np
         v_sum(2) = v_sum(2)/Np
         v_sum(3) = v_sum(3)/Np
         sf = SQRT(init_temp/temp)

         v2_sum = 0
         do ii = 1,Np
           init_vel(ii,:) = (init_vel(ii,:)-v_sum(:))*sf
           v2_sum = v2_sum + init_vel(ii ,1)**2+init_vel(ii ,2)**2+init_vel(ii ,3)**2
         end do
         v2_sum=v2_sum/(3*Np)
         !print*,v2_sum
      end subroutine VELOCITY_ALLOCATION
!***********************************************************************************************!
!SUBROUTINE4:GAUSSIAN DISTRIBUTION                                                              !
!***********************************************************************************************!
      subroutine RAN_GAUSS(ran_gauss1 , ran_gauss2 , ran_gauss3)
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Other_Variables
      implicit none
      integer:: ii                                                    !for test
      real(kind=double)::ran_uniform  , v1 , v2 ,v3 , v4 ,v5 ,v6 ,rsq
      real(kind=double),intent(out)::ran_gauss1 , ran_gauss2 , ran_gauss3
      1  call random_number(v1)
         v1=2*v1 - 1.0
         call random_number(v2)
         v2=2*v2-1.0
         rsq=v1**2+v2**2
         if(rsq .gt. 1.0 .or. rsq .le. 0) GOTO 1
            ran_gauss1 = v1*sqrt(-2.0*log(rsq)/rsq)

      2  call random_number(v3)                                        !Box Muller Algorithm to generate random number gaussian
         v3=2*v3 - 1.0
         call random_number(v4)
         v4=2*v4-1.0
         rsq=v3**2+v4**2
         if(rsq .gt. 1.0 .or. rsq .le. 0) GOTO 2
            ran_gauss2 = v3*sqrt(-2.0*log(rsq)/rsq)

      3  call random_number(v5)
         v5=2*v5 - 1.0
         call random_number(v6)
         v6=2*v6-1.0
         rsq=v5**2+v6**2
         if(rsq .gt. 1.0 .or. rsq .le. 0) GOTO 2
            ran_gauss3 = v5*sqrt(-2.0*log(rsq)/rsq)


      end  subroutine RAN_GAUSS
!*******************************************************************************************!
!SUBROUTINE 4:FORCE CALCULATION                                                             !
!*******************************************************************************************!
      subroutine FORCE
      Use Precision_Control
      Use Parameter_Control
      Use Read_Variables
      Use Other_Variables
      use mpi
      implicit none
      integer ::iparticle , jparticle , iproc
      real(kind=double)::rij(3) ,force_term(3) ,  rsq , r2i , r6i 

!for mpi based parallalization
      integer::an_id
      integer::start_particle , end_particle , deff_to_add ,  start_particle_dummy , end_particle_dummy
      real(kind=double) :: partial_force(Np ,ndim) , local_force(Np , ndim) ,  potential_energy , pe , R_dummy(Np , ndim)
      real(kind=double) ::check_term
      !if(my_id .eq. 0) then 
      !  write(500 , *) R(1 ,1)
      !end if
       num_procs = nop
      if (my_id .eq. root_process) then
       
      
         particle_per_process = int(Np/num_procs)
         particle_for_root    = Np -(particle_per_process*(num_procs-1))
         deff_to_add          = particle_for_root -particle_per_process
         do an_id = 1, num_procs -1
            start_particle = an_id*particle_per_process + 1 +deff_to_add
            end_particle   = start_particle + particle_per_process -1
            if (an_id .eq. num_procs -1) end_particle = Np
!            print*, start_particle , end_particle 
            call MPI_SEND(start_particle , 1 , MPI_INT , an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(end_particle , 1 , MPI_INT , an_id, send_data_tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(R , 3*Np , MPI_DOUBLE_PRECISION , an_id, send_data_tag, MPI_COMM_WORLD, ierr)      
         end do
         net_force = 0
         potential = 0
         
! calculation in root process
         do iparticle = 1 , particle_for_root
            do jparticle = iparticle+1 , Np
               rij(:) = R(iparticle , :) - R(jparticle , :)
               rij(:) = rij(:) -boxlen*anint(rij(:)/boxlen)
               rsq = rij(1)**2+rij(2)**2+rij(3)**2
               if ( rsq .lt. rcut) then
                  r2i = 1/(rsq)
                  r6i = r2i**3
                  if (iparticle .le. Nsolvent .and. jparticle .le. Nsolvent) then
                     epsil = 1.0
                  else if( iparticle .le. Nsolvent .and. jparticle .gt. Nsolvent) then
                         epsil = 0.3
                  else
                         epsil = 0.5
                  end if
                  force_term(:) = 48.0*r2i*epsil*r6i*(r6i-0.50)*rij(:)
                  net_force(iparticle ,:)=net_force(iparticle,:)+force_term(:)
                  net_force(jparticle ,:)=net_force(jparticle,:)-force_term(:)
                  potential = potential+ 4.0*epsil*r6i*(r6i - 1)
               end if
            end do
         end do

!gathering data from processes
        do an_id = 1 , num_procs-1
           call MPI_RECV(partial_force , 3*Np , MPI_DOUBLE_PRECISION , an_id , MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
           call MPI_RECV(pe , 1 , MPI_DOUBLE_PRECISION , an_id , MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
           do iparticle = 1 , Np
           net_force(iparticle , :) = net_force(iparticle , :) +partial_force(iparticle , :)
           end do
           potential = potential +pe
           print* , pe , an_id , potential
        end do
        write(30 ,*) potential
        call MPI_RECV(check_term, 1 , MPI_DOUBLE_PRECISION , 3 , MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      !   write(500 , *) check_term
!work in other nodes
      else
           call MPI_RECV(start_particle_dummy , 1 , MPI_INT , root_process , MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
           call MPI_RECV(end_particle_dummy , 1 , MPI_INT , root_process , MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
           call MPI_RECV(R_dummy , 3*Np , MPI_DOUBLE_PRECISION , root_process , MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
       local_force= 0
       potential_energy = 0
       do iparticle =  start_particle_dummy , end_particle_dummy
          do jparticle = iparticle+1 , Np
             rij(:) = R_dummy(iparticle , :) - R_dummy(jparticle , :)
             rij(:) = rij(:) -boxlen*anint(rij(:)/boxlen)
             rsq = rij(1)**2+rij(2)**2+rij(3)**2
             if ( rsq .lt. rcut) then
                r2i = 1/(rsq)
                r6i = r2i**3
                if (iparticle .le. Nsolvent .and. jparticle .le. Nsolvent) then
                  epsil = 1.0
                else if( iparticle .le. Nsolvent .and. jparticle .gt. Nsolvent) then
                  epsil = 0.3
                else
                  epsil = 0.5
                end if
                force_term(:) = 48.0*r2i*epsil*r6i*(r6i-0.50)*rij(:)
                potential_energy = potential_energy+ 4.0*epsil*r6i*(r6i - 1)
                local_force(iparticle , :) = local_force(iparticle,:)+force_term(:)
                local_force(jparticle,:) = local_force(jparticle,:)-force_term(:)
             end if
          end do
       end do
       do iproc = 1 , num_procs-1
          if(my_id .eq. iproc) then
            call MPI_SEND(local_force , 3*Np , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(potential_energy , 1 , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
          end if
       end do
       !if(my_id .eq. 1) then 
       !  call MPI_SEND(local_force , 3*Np , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
       !  call MPI_SEND(potential_energy , 1 , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
       !else if(my_id .eq. 2) then
       !  call MPI_SEND(local_force , 3*Np , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
       !  call MPI_SEND(potential_energy , 1 , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
       !else if(my_id .eq. 3) then
       !  call MPI_SEND(local_force , 3*Np , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
       !  call MPI_SEND(potential_energy , 1 , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr)
       !end if
!check
       if(my_id .eq. 3) then 
         call MPI_SEND(R_dummy(1:1 , 1) , 1 , MPI_DOUBLE_PRECISION , 0 , return_data_tag, MPI_COMM_WORLD, ierr) 
       end if


      end if
      end subroutine FORCE 
