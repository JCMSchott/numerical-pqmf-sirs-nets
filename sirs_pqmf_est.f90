module pqmf
   use types
   use geraRede

   implicit none
   
   real(dp), allocatable :: I_i(:), S_i(:) ! R(:) = 1 - I(:) - S(:)
   real(dp), allocatable :: RI_ij(:), SI_ij(:), SS_ij(:)

   character(len=500) :: local
   character(len=1000) :: local_arquivo, I_i_arquivo, S_i_arquivo
   character(len=1000) :: SI_ij_arquivo, RI_ij_arquivo, SS_ij_arquivo
   character(len=1000) :: lb_vs_I_i_arquivo, t_vs_I_im_arquivo

   real(dp) :: lambda
   real(dp) :: t
   
   logical :: foi_I_i, foi_SI_ij, foi_RI_ij
   real(dp):: qualvalor
   
   integer  :: nlambda
   real(dp) :: lambda0
   real(dp) :: lambdaf
   real(dp) :: Iest
      
   real(dp) :: dlambda  
   
   contains

   subroutine condicao_inicial(this, flag, file_exists, Nichar, gamachar, lambdachar)     

      class(grafo), intent(in) :: this
      integer:: j1, j12, j21, j2, j23, j3
      integer :: sum_deg
      real(dp) :: sorteio
      real(dp), parameter :: I_0 = 0.01_dp
      logical :: isI_i
      integer :: flag
      logical :: file_exists
      integer :: iost
      character(len=30) :: Nichar
      character(len=10) :: gamachar
      character(len=10) :: lambdachar
      character(len=100) :: arquivo_char2
      
      sum_deg = sum(this%deg)
            
      foi_I_i = .False.
      foi_SI_ij = .False.
      foi_RI_ij = .False.
            
      !######AlocaListas#####################
      if(allocated(I_i)) deallocate(I_i)
      allocate(I_i(this%nodes))
      I_i = 0.0_dp

      if(allocated(S_i)) deallocate(S_i)
      allocate(S_i(this%nodes))
      S_i = 0.0_dp

      if(allocated(RI_ij)) deallocate(RI_ij)
      allocate(RI_ij(sum_deg))      
      RI_ij = 0.0_dp

      if(allocated(SI_ij)) deallocate(SI_ij)
      allocate(SI_ij(sum_deg))      
      SI_ij = 0.0_dp

      if(allocated(SS_ij)) deallocate(SS_ij)
      allocate(SS_ij(sum_deg))
      SS_ij = 0.0_dp
      
      !######AlocaListas#####################
      t = 0.0_dp
      
      write(*,*) "Flag= ", flag
      if(flag == 1)then
         
!        open(unit=26, file=trim(adjustl(I_i_arquivo)), status='unknown')
         
         if(file_exists)then
            close(20)
            open(unit = 20, file = trim(adjustl(lb_vs_I_i_arquivo)), status = 'unknown')

            write(*,*) "Arquivos em lbd ja existem"
            !t = 0.0_dp


le2:        do
               read(20,*, iostat = iost) lambda, Iest               
               if(iost > 0)then
                  write(*,*) "Erro incompreendido"
                  exit le2
               elseif(iost < 0)then
                  write(*,*) "Fim do arquivo alcancado"
                  exit le2
               endif
               write(*,*) "Valor atual de lambda ", lambda
            enddo le2
!           lambda = lambda + dlambda
            close(20)
            open(unit = 20, file = trim(adjustl(lb_vs_I_i_arquivo)), access = 'append', status = 'unknown')
            write(*,*) "Valor de lambda ", lambda

!############################################################################################
            write(lambdachar, '(f7.4)') lambda
            arquivo_char2 = trim(adjustl('t_vs_I_im_tam_'//trim(adjustl(Nichar))//'_g_'//trim(adjustl(gamachar))//'_lbd_'//trim(adjustl(lambdachar))))
            
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char2))//'.dat'
            t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
            open(unit = 31, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')!############################################################################################            
         else
            stop "Arquivo lbd_vs_I_im... nao existe."
         endif



le1:     do
            read(31,*, iostat = iost) t, Iest
            if(iost > 0)then
               write(*,*) "Erro incompreendido na leitura do tempo"
               exit le1
            elseif(iost < 0)then
               write(*,*) "Fim do arquivo alcancado"
               exit le1
            endif
            write(*,*) "Valor atual de t eh ", t
         enddo le1
         rewind(31)

                  
         do j1 = 1, size(I_i)
            read(26,*) I_i(j1)
         enddo

         rewind(26)
         
         open(unit=27, file=trim(adjustl(S_i_arquivo)), status='unknown')

         do j1 = 1, size(S_i)
            read(27,*) S_i(j1)
         enddo
         
         rewind(27)

         open(unit=28, file=trim(adjustl(SI_ij_arquivo)), status='unknown')

         do j1 = 1, size(SI_ij)
            read(28,*) SI_ij(j1)
         enddo
         
         rewind(28)

         open(unit=29, file=trim(adjustl(RI_ij_arquivo)), status='unknown')

         do j1 = 1, size(RI_ij)
            read(29,*) RI_ij(j1)
         enddo
         
         rewind(29)

         open(unit=30, file=trim(adjustl(SS_ij_arquivo)), status='unknown')

         do j1 = 1, size(SS_ij)
            read(30,*) SS_ij(j1)
         enddo
         rewind(30)
         
      elseif(flag == 0)then
         if(file_exists)then
            close(20)
            open(unit = 20, file = trim(adjustl(lb_vs_I_i_arquivo)), status = 'unknown')

            write(*,*) "Arquivos ja existem"
            t = 0.0_dp
            
le3:        do
               read(20,*, iostat = iost) lambda, Iest               
               if(iost > 0)then
                  write(*,*) "Erro incompreendido"
                  exit le3
               elseif(iost < 0)then
                  write(*,*) "Fim do arquivo alcancado"
                  exit le3
               endif
               write(*,*) "Valor atual de lambda ", lambda
            enddo le3
            lambda = lambda + dlambda
            close(20)
            open(unit = 20, file = trim(adjustl(lb_vs_I_i_arquivo)), access = 'append', status = 'unknown')
            write(*,*) "Valor de lambda ", lambda
         else
            write(*,*) "Criando arquivos"     
         endif
               
         do j1 = 1, this%nodes
            I_i(j1) = I_0
            S_i(j1) = 1.0_dp - I_i(j1)
         enddo

         do j1 = 1, this%nodes
                  
            do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
               if(SI_ij(j12) > 0.0_dp) cycle
               j2 = this%listAdj(j12)
             
               do j23 = this%aux(j2), this%aux(j2) + this%deg(j2) - 1
                  if(this%listAdj(j23) == j1)then
                     j21 = j23
                     exit
                  endif
               enddo
                                       
               SI_ij(j12) = S_i(j1) * I_i(j2)
               
               SI_ij(j21) = S_i(j2) * I_i(j1)
 
               SS_ij(j12) = S_i(j1) * S_i(j2)
            
               SS_ij(j21) = SS_ij(j12)
                                 
            enddo
         enddo
      else
         stop "Nao foi possivel settar as condicoes iniciais com o valor indicado de flag."   
      endif
                         
   end subroutine       
     
   subroutine k4_sirs_pqmf(dt, t, this, ms, ns, alp, lbd, mu, I_i, S_i, RI_ij, SI_ij, SS_ij, fI_i, fS_i, fRI_ij, fSI_ij, fSS_ij)
      
      use types
      use geraRede
      
      real(dp) :: dt
      real(dp) :: t
      class(grafo), intent(in) :: this         
      integer :: ms
      integer :: ns
      real(dp) :: alp, lbd, mu

      real(dp) :: I_i(ms)
      real(dp) :: S_i(ms)
      real(dp) :: RI_ij(ns)
      real(dp) :: SI_ij(ns)
      real(dp) :: SS_ij(ns)


      real(dp) :: I_ip(ms)
      real(dp) :: S_ip(ms)
      
      real(dp) :: k1_I_i(ms), k2_I_i(ms), k3_I_i(ms), k4_I_i(ms)
      real(dp) :: k1_S_i(ms), k2_S_i(ms), k3_S_i(ms), k4_S_i(ms)

      real(dp) :: RI_ijp(ns)
      real(dp) :: SI_ijp(ns)
      real(dp) :: SS_ijp(ns)
      
      real(dp) :: k1_RI_ij(ns), k2_RI_ij(ns), k3_RI_ij(ns), k4_RI_ij(ns)
      real(dp) :: k1_SI_ij(ns), k2_SI_ij(ns), k3_SI_ij(ns), k4_SI_ij(ns)
      real(dp) :: k1_SS_ij(ns), k2_SS_ij(ns), k3_SS_ij(ns), k4_SS_ij(ns)

      real(dp) :: argI_i(ms), argS_i(ms), argRI_ij(ns), argSI_ij(ns), argSS_ij(ns)

      integer :: l1

      !#########################INTERFACE#############################
      
      interface
            
 	function fI_i(t, this, ms, ns, lbd, mu, I_i, SI_ij)
            use types
            use geraRede
            
            real(dp), intent(in) :: t
            class(grafo), intent(in) :: this       
            integer, intent(in) :: ms
            integer, intent(in) :: ns
            real(dp), intent(in) :: lbd, mu
            real(dp) :: I_i(ms)
            real(dp) :: SI_ij(ns)
            real(dp) :: soma       
            real(dp) :: fI_i(ms)
            integer :: k12 
 	end function      


 	function fS_i(t, this, ms, ns, alp, lbd, I_i, S_i, SI_ij)
            use types
            use geraRede
            
            real(dp), intent(in) :: t
            class(grafo), intent(in) :: this         
            integer, intent(in) :: ms
            integer, intent(in) :: ns
            real(dp), intent(in) :: alp, lbd
            real(dp) :: I_i(ms)
            real(dp) :: S_i(ms)
            real(dp) :: SI_ij(ns)       
            real(dp) :: fS_i(ms)
            real(dp) :: soma
            integer :: k12

 	end function 			
 
        function fRI_ij(t, this, ms, ns, alp, lbd, mu, I_i, S_i, RI_ij, SI_ij, SS_ij)
            use types
            use geraRede

            real(dp), intent(in) :: t
            class(grafo), intent(in) :: this         
            integer, intent(in) :: ms
            integer, intent(in) :: ns
            real(dp), intent(in) :: alp, lbd, mu
            real(dp) :: I_i(ms)
            real(dp) :: S_i(ms)
            real(dp) :: RI_ij(ns)
            real(dp) :: SI_ij(ns)
            real(dp) :: SS_ij(ns)           
            real(dp) :: fRI_ij(ns)

            real(dp) :: dum_II_ij, dum_RS_ijSj, dum_SS_ijSj
       
            integer :: edge_ji
            integer :: viz1
            integer :: k3    ! k1, k2 ja existem no geraRede
            integer :: k12, k12_dum, k21, k23, k32  
            real(dp) :: soma

        end function  

         function fSI_ij(t, this, ms, ns, alp, lbd, mu, S_i, RI_ij, SI_ij, SS_ij)

            use types
            use geraRede

            real(dp), intent(in) :: t
            class(grafo), intent(in) :: this         
            integer, intent(in) :: ms
            integer, intent(in) :: ns
            real(dp), intent(in) :: alp, lbd, mu
            real(dp) :: S_i(ms)
            real(dp) :: RI_ij(ns)
            real(dp) :: SI_ij(ns)
            real(dp) :: SS_ij(ns)           
            real(dp) :: fSI_ij(ns) 			
           
            real(dp) :: dum_SI_ijSi, dum_SS_ijSj
       
            integer :: k3 
            integer :: k12, k12_dum, k21, k23, k32  
            real(dp) :: soma

 	end function


        function fSS_ij(t, this, ms, ns, alp, lbd, mu, S_i, SI_ij, SS_ij)

            use types
            use geraRede

            real(dp), intent(in) :: t
            class(grafo), intent(in) :: this         
            integer, intent(in) :: ms
            integer, intent(in) :: ns
            real(dp), intent(in) :: alp, lbd, mu
            real(dp) :: S_i(ms)
            real(dp) :: SI_ij(ns)
            real(dp) :: SS_ij(ns)           
            real(dp) :: fSS_ij(ns)
            real(dp) :: dum_II_ij, dum_RS_ij, dum_RS_ji, dum_SS_ijSi, dum_SS_ijSj
       
            integer :: edge_ji
            integer :: viz1
            integer :: k3   !k1, k2 ja existem
            integer :: k12, k12_dum, k21, k23, k32  
            real(dp) :: soma

 	end function 			      
      end interface        
      !#########################INTERFACE#############################
      
      !#########################SUBROTINA#############################

!#################################################################################################
		k1_I_i = dt * fI_i(t, this, ms, ns, lbd, mu, I_i, SI_ij)		
		
		k1_S_i = dt * fS_i(t, this, ms, ns, alp, lbd, I_i, S_i, SI_ij)
		
		k1_RI_ij = dt * fRI_ij(t, this, ms, ns, alp, lbd, mu, I_i, S_i, RI_ij, SI_ij, SS_ij)
		
		k1_SI_ij = dt * fSI_ij(t, this, ms, ns, alp, lbd, mu, S_i, RI_ij, SI_ij, SS_ij)
		
		k1_SS_ij = dt * fSS_ij(t, this, ms, ns, alp, lbd, mu, S_i, SI_ij, SS_ij)
!#################################################################################################
		!###############################		
                argI_i = I_i + 0.5_dp * k1_I_i
                argS_i = S_i + 0.5_dp * k1_S_i
                argRI_ij = RI_ij + 0.5_dp * k1_RI_ij
                argSI_ij = SI_ij + 0.5_dp * k1_SI_ij			
                argSS_ij = SS_ij + 0.5_dp * k1_SS_ij
		!###############################
!#################################################################################################
		k2_I_i = dt * fI_i(t + (dt/2_dp), this, ms, ns, lbd, mu, argI_i, argSI_ij)		
				
		k2_S_i = dt * fS_i(t + (dt/2_dp), this, ms, ns, alp, lbd, argI_i, argS_i, argSI_ij)

		k2_RI_ij = dt * fRI_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k2_SI_ij = dt * fSI_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
		                
		k2_SS_ij = dt * fSS_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argS_i, argSI_ij, argSS_ij)
!#################################################################################################
		!###############################
		argI_i = I_i + 0.5_dp * k2_I_i
		argS_i = S_i + 0.5_dp * k2_S_i
		argRI_ij = RI_ij + 0.5_dp * k2_RI_ij
                argSI_ij = SI_ij + 0.5_dp * k2_SI_ij			
		argSS_ij = SS_ij + 0.5_dp * k2_SS_ij
		!###############################
!#################################################################################################
		k3_I_i = dt * fI_i(t + (dt/2_dp), this, ms, ns, lbd, mu, argI_i, argSI_ij)		
		
		k3_S_i = dt * fS_i(t + (dt/2_dp), this, ms, ns, alp, lbd, argI_i, argS_i, argSI_ij)
		
		k3_RI_ij = dt * fRI_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
		                
		k3_SI_ij = dt * fSI_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k3_SS_ij = dt * fSS_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argS_i, argSI_ij, argSS_ij)
!#################################################################################################
               !###############################		
               argI_i = I_i + k3_I_i
               argS_i = S_i + k3_S_i
               argRI_ij = RI_ij + k3_RI_ij
               argSI_ij = SI_ij + k3_SI_ij			
               argSS_ij = SS_ij + k3_SS_ij
               !###############################
!#################################################################################################
		k4_I_i = dt * fI_i(t + dt, this, ms, ns, lbd, mu, argI_i, argSI_ij)		
		
		k4_S_i = dt * fS_i(t + dt, this, ms, ns, alp, lbd, argI_i, argS_i, argSI_ij)
		
		k4_RI_ij = dt * fRI_ij(t + dt, this, ms, ns, alp, lbd, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k4_SI_ij = dt * fSI_ij(t + dt, this, ms, ns, alp, lbd, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k4_SS_ij = dt * fSS_ij(t + dt, this, ms, ns, alp, lbd, mu, argS_i, argSI_ij, argSS_ij)
!################################################################################################
!               do l1 = 1, this%nodes
!                 I_i(l1) = I_i(l1) + 1d0/6d0 * (k1_I_i(l1) + 2d0 * k2_I_i(l1) + 2d0 * k3_I_i(l1) + k4_I_i(l1))
!                  if(I_i(l1) < 0.0_dp)then
!                     foi_I_i = .True.
!                     qualvalor = I_i(l1)
!                  endif
!               enddo

!               do l1 = 1, sum(this%deg)
!                  RI_ij(l1) = RI_ij(l1) + 1d0/6d0 * (k1_RI_ij(l1) + 2d0 * k2_RI_ij(l1) + 2d0 * k3_RI_ij(l1) + k4_RI_ij(l1))
                  
!                  if(RI_ij(l1) < 0.0_dp)then
!                     foi_RI_ij = .True.
!                     qualvalor = RI_ij(l1)
!                  endif
!               enddo
               
!               do l1 = 1, sum(this%deg)
!                  SI_ij(l1) = SI_ij(l1) + 1d0/6d0 * (k1_SI_ij(l1) + 2d0 * k2_SI_ij(l1) + 2d0 * k3_SI_ij(l1) + k4_SI_ij(l1))
                  
!                  if(SI_ij(l1) < 0.0_dp)then
!                     foi_SI_ij = .True.
!                     qualvalor = SI_ij(l1)
!                  endif
!               enddo
!################################################################################################
               I_i = I_i + 1d0/6d0 * (k1_I_i + 2d0 * k2_I_i + 2d0 * k3_I_i + k4_I_i)

               S_i = S_i + 1d0/6d0 * (k1_S_i + 2d0 * k2_S_i + 2d0 * k3_S_i + k4_S_i)

               RI_ij = RI_ij + 1d0/6d0 * (k1_RI_ij + 2d0 * k2_RI_ij + 2d0 * k3_RI_ij + k4_RI_ij)

               SI_ij = SI_ij + 1d0/6d0 * (k1_SI_ij + 2d0 * k2_SI_ij + 2d0 * k3_SI_ij + k4_SI_ij)
               
               SS_ij = SS_ij + 1d0/6d0 * (k1_SS_ij + 2d0 * k2_SS_ij + 2d0 * k3_SS_ij + k4_SS_ij)
!#########################SUBROTINA##############################################################
                                             
   end subroutine
   
end module

!#######################


program sirs_pqmp_est

   !#####preambulo######   

!  use mod_rndgen
   use geraRede
   use mod_tools
   use types
   use pqmf
   !#####preambulo######

   implicit none

!  real(dp), allocatable :: I_i(:), S_i(:) ! R(:) = 1 - I(:) - S(:)
!  real(dp), allocatable :: RI_ij(:), SI_ij(:), SS_ij(:)
      
   integer, parameter :: t0 = 0.0_dp
   real(dp), parameter :: tf = 1000.0_dp
   real(dp), parameter :: dt = 0.01_dp
   integer, parameter :: npts = int((tf - t0)/dt)
   
   
   real(dp), parameter :: alp = 0.5_dp
   real(dp), parameter :: mu = 1.0_dp
   integer :: semente
   
!  type(rndgen) :: gerador
   
   real(dp), allocatable :: exp_gama(:)
   integer :: i1, i2, i3, i4, i5, i6
   integer :: j1
   integer :: sum_deg
   integer :: k_i
   real(dp) :: k_f
   
   type(grafo_PL_UCM) :: rede
   
   real :: I_im0
   real(dp) :: R_im0
   real(dp) :: RI_ijm0
   real(dp) :: S_im0
   real(dp) :: SI_ijm0
   real(dp) :: SS_ijm0
   
   real :: I_im
   real(dp) :: R_im
   real(dp) :: RI_ijm
   real(dp) :: S_im
   real(dp) :: SI_ijm
   real(dp) :: SS_ijm
   
   integer :: I_im_repetido
   real(dp) :: delta_I_im
   real(dp) :: toll
   logical :: tem_trans
   
   character(len=100) :: arquivo_char1
   real(dp), parameter :: t_rec = 1.0_dp * int(tf/500.0_dp)
   real(dp) :: t_rec1
   logical :: arquivo_existe
   integer :: n_args
   integer :: flag1
   character(len=1) :: flag_char
   character(len=10) :: lambda_char 
   integer, allocatable :: Ni(:)
   character(len=30), allocatable :: Ni_char(:)
   character(len=10) :: gama_char
   
   dlambda = 1.0_dp * (lambdaf - lambda0)/nlambda
   
!   local = trim(adjustl('/home/jota/S_comp/S_num/SIRS_pQMF/Rst/'))
   local = trim(adjustl('/home/jota/SIRS_pQMF/Rst/'))   
!##########################################################################################
!		Aloca lista de tamanhos
!##########################################################################################
   allocate(Ni(8))

   allocate(Ni_char(8))

   allocate(exp_gama(3))
	
   exp_gama(1) = 2.3_dp; exp_gama(2) = 2.7_dp; exp_gama(3) = 3.5_dp
	
   Ni(1) = 1000; Ni(2) = 3000; Ni(3) = 10000; Ni(4) = 30000
   Ni(5) = 100000; Ni(6) = 300000; Ni(7) = 1000000; Ni(8) = 3000000
!  Ni(9) = 10000000; Ni(10) = 30000000

   Ni_char(1) = '1k'; Ni_char(2) = '3k'; Ni_char(3) = '10k'; Ni_char(4) = '30k'
   Ni_char(5) = '100k'; Ni_char(6) = '300k'; Ni_char(7) = '1M'; Ni_char(8) = '3M'
!  Ni_char(9) = '10M'; Ni_char(10) = '30M'

!################################################################################   
! Defino os parametros lambda.
!  12_09_2020:   01:21am. refinamento.   
   nlambda = 1000
   
   lambda0 = 5.0d-3 
   lambdaf = 1.0_dp
   dlambda = 1.0_dp * (lambdaf - lambda0)/nlambda
!################################################################################   
   semente = 995887671
   
   k_i = 3
!################################################################################   
   n_args = iargc()

   if(n_args == 1)then
      call getarg(n_args, flag_char)
      read(flag_char,*) flag1
      if((flag1 /= 1).and.(flag1 /=0)) stop "Argumento invalido"
   else
      stop "Numero de argumentos invalido"
   endif   

!  call gerador%init(semente)
   
   do i1 = 1, 1 !size(exp_gama)
      do i2 = 1, 1 ! size(Ni)
         
         toll = 1.0_dp/Ni(i2)
         
         if(exp_gama(i1) > 3.0_dp)then
            k_f = 1.0_dp * (1.0_dp * Ni(i2))**(1.0_dp/(exp_gama(i1) -1.0_dp))
         else
            k_f = 1.0_dp * sqrt(1.0_dp * Ni(i2))
         endif         
         
         call rede%iniciaGrafo(Ni(i2))
         call rede%inicia(k_i, k_f, exp_gama(i1), semente)   
         call rede%liga(semente, .False.)
         
         sum_deg = sum(rede%deg)
         !############CondicaoInicial#################
            
         !############CondicaoInicial#################
         
         !############Dinamica#################
         
!        tem_trans = .False.

         lambda = lambda0
  
         Ni_char(i2) = trim(adjustl(Ni_char(i2)))
            
!        write(gama_char, '(f4.2)') exp_gama(i1)
         write(gama_char, '(f5.2)') exp_gama(i1)         
         
         gama_char = trim(adjustl(gama_char))


!############################################################################################
         arquivo_char1 = trim(adjustl('t_vs_I_i_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
         
         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
                  
         I_i_arquivo = local_arquivo
         
         open(unit = 26, file = trim(adjustl(I_i_arquivo)), status ='unknown')
                     
!############################################################################################
         arquivo_char1 = trim(adjustl('t_vs_S_i_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))

         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            
         S_i_arquivo = local_arquivo
         open(unit = 27, file = trim(adjustl(S_i_arquivo)), status ='unknown')
!############################################################################################
            
         arquivo_char1 = trim(adjustl('t_vs_SI_ij_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
         
         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            
         SI_ij_arquivo = local_arquivo
            
         open(unit = 28, file = trim(adjustl(SI_ij_arquivo)), status = 'unknown')
         
!############################################################################################
         arquivo_char1 = trim(adjustl('t_vs_RI_ij_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
         RI_ij_arquivo = local_arquivo
            
         open(unit = 29, file = trim(adjustl(RI_ij_arquivo)), status = 'unknown')
         
!############################################################################################
         arquivo_char1 = trim(adjustl('t_vs_SS_ij_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
            
         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            
         SS_ij_arquivo = local_arquivo
            
         open(unit = 30, file = trim(adjustl(SS_ij_arquivo)), status = 'unknown')

!        endif
            
!############################################################################################
         arquivo_char1 = trim(adjustl('lbd_vs_I_im_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))

         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
         lb_vs_I_i_arquivo = local_arquivo
         
         lb_vs_I_i_arquivo = trim(adjustl(lb_vs_I_i_arquivo))
         
         inquire(file = trim(adjustl(lb_vs_I_i_arquivo)), exist = arquivo_existe)
         
!        if(.not.arquivo_existe)then
         
         open(unit = 20, file = trim(adjustl(lb_vs_I_i_arquivo)), access = 'append', status = 'unknown')

!############################################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_S_im_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))

            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
            open(unit = 21, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
!############################################################################################

!############################################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_R_im_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
         
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
            open(unit = 22, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')

!############################################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_SI_ijm_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
         
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat' 
         
            open(unit = 23, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')

!############################################################################################

            arquivo_char1 = trim(adjustl('lbd_vs_RI_ijm_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))
         
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'

            open(unit = 24, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
!############################################################################################

            arquivo_char1 = trim(adjustl('lbd_vs_SS_ijm_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))

            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
            open(unit = 25, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')

            arquivo_char1 = trim(adjustl('lbd_vs_t_conv_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))))

            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'


            open(unit = 32, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
         
!        endif
!#############################################################################################
         
         do i3 = 1, nlambda
            
            call condicao_inicial(rede, flag1, arquivo_existe, Ni_char(i2), gama_char, lambda_char) 
                                    
            if(flag1 == 0)then
!############################################################################################
               write(lambda_char, '(f7.4)') lambda
               arquivo_char1 = trim(adjustl('t_vs_I_im_tam_'//trim(adjustl(Ni_char(i2)))//'_g_'//trim(adjustl(gama_char))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 31, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
!############################################################################################
            endif
            !############TempoRoda#################

                        
            I_im = 1.0*(1.0_dp * sum(I_i)/size(I_i))
            S_im = 1.0_dp * sum(S_i)/size(S_i)
            R_im = 1.0_dp - I_im - S_im
            SI_ijm = 1.0_dp * sum(SI_ij)/size(SI_ij)
            SS_ijm = 1.0_dp * sum(SS_ij)/size(SS_ij)
            RI_ijm = 1.0_dp * sum(SI_ij)/size(RI_ij)            
           
            I_im_repetido = 0
            
            t_rec1 = t_rec
                     
din1:       do i4 = 1, npts
                           
               I_im0 = 1.0*I_im                              
               S_im0 = S_im
               R_im0 = 1.0_dp - S_im0 - I_im0
               
               RI_ijm0 = RI_ijm
               SI_ijm0 = SI_ijm
               SS_ijm0 = SS_ijm
               
               
               write(31,*) t, I_im
                              
               call k4_sirs_pqmf(dt, t, rede, rede%nodes, sum_deg, alp, lambda, mu, I_i, S_i, RI_ij, SI_ij, SS_ij, fI_i, fS_i, fRI_ij, fSI_ij, fSS_ij)
               
                              
               I_im = 1.0*(1.0_dp * sum(I_i)/size(I_i))
               S_im = 1.0_dp * sum(S_i)/size(S_i)
               
               R_im = 1.0_dp - I_im - S_im
               
               RI_ijm = 1.0_dp * sum(RI_ij)/size(RI_ij)                           
               SI_ijm = 1.0_dp * sum(SI_ij)/size(SI_ij)
               SS_ijm = 1.0_dp * sum(SS_ij)/size(SS_ij)
                                              
               !################SeDaNegativoOuRepete##################

!              if((foi_I_i).or.(foi_SI_ij).or.(foi_RI_ij))then 
            
!                  I_im = I_im0               
!                  R_im = R_im0
!                  S_im = S_im0
!                  RI_ijm = RI_ijm0
!                  SI_ijm = SI_ijm0
!                  SS_ijm = SS_ijm0
                  
!                  if(foi_I_i)then
!                     write(*,*) "I_i deu valor negativo com o valor ", qualvalor
!                  elseif(foi_RI_ij)then
!                     write(*,*) "RI_ij deu valor negativo com o valor ", qualvalor
!                  elseif(foi_SI_ij)then
!                     write(*,*) "SI_ij deu valor negativo com o valor ", qualvalor
!                  endif                  
!                  write(*,*) "No instante de tempo t ", t   
!                  exit din1           
                 if(I_im == I_im0)then
                  if(I_im == I_im0) I_im_repetido = I_im_repetido + 1
                  if(I_im_repetido == 5)then
                     I_im = I_im0               
                     R_im = R_im0
                     S_im = S_im0
                     RI_ijm = RI_ijm0
                     SI_ijm = SI_ijm0
                     SS_ijm = SS_ijm0
                     write(*,*) "Repetiu o mesmo valor varias vezes"                     
                     exit din1
                  endif 
               endif
               !################SeDaNegativoOuRepete##################
               
               if(I_im <= 1.0_dp/(5.0_dp * rede%nodes))then
                  I_im = 0.0_dp
                  S_im = 1.0_dp
                  R_im = 0.0_dp
                  RI_ijm = 0.0_dp
                  SI_ijm = 0.0_dp
                  SS_ijm = 1.0_dp
                  exit din1
               endif
                              
               t = t + dt

               if( (t < t_rec1 + 1.5_dp * dt) .and. (t > t_rec1 - 1.5_dp * dt))then
            
                   write(26,*) t
                   
                   write(26,*) lambda
                          
                   do j1 = 1, size(I_i)
                      write(26,*) I_i(j1)
                   enddo
                   rewind(26)
                   
                   do j1 = 1, size(S_i)
                      write(27,*) S_i(j1)
                   enddo
                   rewind(27)
                   
                   do j1 = 1, size(SI_ij)
                      write(28,*) SI_ij(j1)
                   enddo
                   rewind(28)
                   
                   do j1 = 1, size(RI_ij)
                      write(29,*) RI_ij(j1)
                   enddo
                   rewind(29)
                   
                   do j1 = 1, size(SS_ij)
                      write(30,*) SS_ij(j1)
                   enddo
                   rewind(30)
                            
                   t_rec1 = t_rec1 + t_rec
               endif
                
            enddo din1

            write(20, *) lambda, I_im
            write(21, *) lambda, S_im
            write(22, *) lambda, R_im
            write(23, *) lambda, SI_ijm
            write(24, *) lambda, RI_ijm
            write(25, *) lambda, SS_ijm
            write(32, *) lambda, t
            
            lambda = lambda + dlambda
            
         enddo
         !############Dinamica#################
         do i5 = 20, 31
            close(i5)
         enddo
         
      enddo
   enddo
   
!  stop
   
   contains

   function fI_i(t, this, m_sitios, n_arestas, lambda, mu, I_i, SI_ij)

      use types
      use geraRede

      real(dp), intent(in) :: t
      class(grafo), intent(in) :: this       
      integer, intent(in) :: m_sitios
      integer, intent(in) :: n_arestas
      real(dp), intent(in) :: lambda, mu
       
      real(dp) :: I_i(m_sitios)
      real(dp) :: SI_ij(n_arestas)
       
      real(dp) :: soma
      real(dp) :: fI_i(m_sitios)
       
      integer :: k12


      do k1 = 1, this%nodes
         fI_i(k1) = -1.0_dp * mu * I_i(k1)
          
          
         soma = 0.0_dp
          
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            soma = soma + SI_ij(k12)
         enddo
         fI_i(k1) = fI_i(k1) + lambda * soma
      enddo
       
!     return
       
   end function
    
       
   function fS_i(t, this, m_sitios, n_arestas, alp, lambda, I_i, S_i, SI_ij)

      use types
      use geraRede

      real(dp), intent(in) :: t
      class(grafo), intent(in) :: this         
      integer, intent(in) :: m_sitios
      integer, intent(in) :: n_arestas
      real(dp), intent(in) :: alp, lambda
      real(dp) :: I_i(m_sitios)
      real(dp) :: S_i(m_sitios)
      real(dp) :: SI_ij(n_arestas)       
      real(dp) :: fS_i(m_sitios)
      real(dp) :: soma
!     integer :: k1
      integer :: k12
       
              
      do k1 = 1, this%nodes 

         fS_i(k1) = alp * (1.0_dp -S_i(k1) -I_i(k1))

         soma = 0.0_dp
          
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            soma = soma + SI_ij(k12)
         enddo
         fS_i(k1) = fS_i(k1) - lambda * soma    
      enddo
       
!     return
       
   end function

   function fRI_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, I_i, S_i, RI_ij, SI_ij, SS_ij)

      use types
      use geraRede

      real(dp), intent(in) :: t
      class(grafo), intent(in) :: this         
      integer, intent(in) :: m_sitios
      integer, intent(in) :: n_arestas
      real(dp), intent(in) :: alp, lambda, mu
      real(dp) :: I_i(m_sitios)
      real(dp) :: S_i(m_sitios)
      real(dp) :: RI_ij(n_arestas)
      real(dp) :: SI_ij(n_arestas)
      real(dp) :: SS_ij(n_arestas)           
      real(dp) :: fRI_ij(n_arestas)
           
      real(dp) :: dum_II_ij, dum_RS_ijSj, dum_SS_ijSj
       
      integer :: edge_ji
      integer :: viz1
      integer :: k3    ! k1, k2 ja existem no geraRede
      integer :: k12, k12_dum, k21, k23, k32  
      real(dp) :: soma

      do k1 = 1, this%nodes
        do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
           k2 = this%listAdj(k12)
           
           do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
              if(this%listAdj(k23) == k1)then
                 k21 = k23
                 exit
              endif
           enddo

!          dum_II_ij = I_i(k1) - SI_ij(k21) - RI_ij(k21)
           
           dum_II_ij = I_i(k2) - SI_ij(k12) - RI_ij(k12)
           
           fRI_ij(k12) = -1.0_dp * (alp + mu) * RI_ij(k12) + mu * dum_II_ij
                      
           soma = 0.0_dp
           do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
              if(k23 == k21) cycle
              soma = soma + SI_ij(k23)
           enddo
           
!          dum_RS_ijSj = (S_i(k2) -SI_ij(k21) - SS_ij(k12))/S_i(k2)   !anterior
           dum_RS_ijSj = (S_i(k2) -SI_ij(k21) - SS_ij(k21))/S_i(k2)                      
           fRI_ij(k12) = fRI_ij(k12) + lambda * dum_RS_ijSj * soma
        enddo
      enddo
      
      return
      
   end function
    
   function fSI_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, S_i, RI_ij, SI_ij, SS_ij)

      use types
      use geraRede

      real(dp), intent(in) :: t
      class(grafo), intent(in) :: this         
      integer, intent(in) :: m_sitios
      integer, intent(in) :: n_arestas
      real(dp), intent(in) :: alp, lambda, mu
      real(dp) :: S_i(m_sitios)
      real(dp) :: RI_ij(n_arestas)
      real(dp) :: SI_ij(n_arestas)
      real(dp) :: SS_ij(n_arestas)           
      real(dp) :: fSI_ij(n_arestas)
           
      real(dp) :: dum_SI_ijSi, dum_SS_ijSj ! dum_II_ij, dum_RS_ij,
       
      integer :: k3     !k1, k2 ja existem 
      integer :: k12, k12_dum, k21, k23, k32  
      real(dp) :: soma
      
      do k1 = 1, this%nodes       
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1 !a1
            k2 = this%listAdj(k12)            

            fSI_ij(k12) = -(mu + lambda) * SI_ij(k12) + alp * RI_ij(k12) 
            
            
            soma = 0.0_dp
            
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1 !k3
               if(this%listAdj(k23) == k1)cycle
               soma = soma + SI_ij(k23)
            enddo

            dum_SS_ijSj = SS_ij(k12)/S_i(k2)
            
            fSI_ij(k12) = fSI_ij(k12) + lambda * dum_SS_ijSj * soma            
            
            
            soma = 0.0_dp
            
            do k12_dum = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
               if(k12_dum == k12) cycle
               soma = soma + SI_ij(k12_dum)
            enddo

            dum_SI_ijSi = SI_ij(k12)/S_i(k1)
            
            fSI_ij(k12) = fSI_ij(k12) - lambda * dum_SI_ijSi * soma
         enddo
         
      enddo
      
      return
      
   end function

   function fSS_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, S_i, SI_ij, SS_ij)

      use types
      use geraRede

      real(dp), intent(in) :: t
      class(grafo), intent(in) :: this         
      integer, intent(in) :: m_sitios
      integer, intent(in) :: n_arestas
      real(dp), intent(in) :: alp, lambda, mu
      real(dp) :: S_i(m_sitios)
      real(dp) :: SI_ij(n_arestas)
      real(dp) :: SS_ij(n_arestas)           
      real(dp) :: fSS_ij(n_arestas)
           
      real(dp) :: dum_II_ij, dum_RS_ij, dum_RS_ji, dum_SS_ijSi, dum_SS_ijSj
       
      integer :: edge_ji
      integer :: viz1
      integer :: k3   !k1, k2 ja existem
      integer :: k12, k12_dum, k21, k23, k32  
      real(dp) :: soma

      do k1 = 1, this%nodes
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) -1
            k2 = this%listAdj(k12)
            
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) -1
               if(this%listAdj(k23) == k1)then
                  k21 = k23
                  exit
               endif
            enddo
                        
            dum_RS_ij = S_i(k2) - SI_ij(k21) - SS_ij(k12)
            dum_RS_ji = S_i(k1) - SI_ij(k12) - SS_ij(k21)
            
            fSS_ij(k12) = alp * (dum_RS_ij + dum_RS_ji)
            
            
            soma = 0.0_dp
            do k12_dum = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
               if(k12_dum == k12) cycle
               soma = soma + SI_ij(k12_dum)
            enddo

            dum_SS_ijSi = SS_ij(k12)/S_i(k1)

            fSS_ij(k12) = fSS_ij(k12) - lambda * dum_SS_ijSi * soma
            
            
            soma = 0.0_dp
            
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               if(k23 == k21) cycle
               soma = soma + SI_ij(k23)
            enddo

            dum_SS_ijSj = SS_ij(k12)/S_i(k2)

            fSS_ij(k12) = fSS_ij(k12) - lambda * dum_SS_ijSj * soma             
         enddo
      
      enddo
      
      return
      
   end function
   
end program
