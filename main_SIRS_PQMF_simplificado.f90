module limiaresCampoMedio
   
   use geraRede
   use mod_rndgen
   
   implicit none
   
   integer, private, parameter :: dp = kind(0.0d0)
   
   real(dp), allocatable :: x(:), y(:)
   real(dp), allocatable :: fi(:)
   real(dp) :: Y4
   
   logical :: salvouCentralidade
   logical :: aberto
   contains
   
   subroutine metodo_bissecao(this, alp, mu, lambda_mean, tole)
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp
      real(dp), intent(in) :: mu
      real(dp), intent(inout) :: lambda_mean
      real(dp) :: lambda_mean_0
      real(dp), intent(in) :: tole
      real(dp) :: lambda_minus, lambda_plus
            
      real(dp) :: dlambda, dlambda0
      real(dp) :: lambda0, lambda_backup

      real(dp) :: autovalor, autovalor0
      real(dp) :: autovalor_minus, autovalor_plus
      real(dp) :: x_norm, y_norm
      integer :: i1
      
      
      lambda0 = limiarQMF(this)
   
      dlambda0 = lambda0/2.0d0
      
      dlambda = dlambda0
     
      !lambda_backup = lambda0
      open(789, file = 'nit_vs_lbd_bissec.dat')
   !====================================================================   
   i1 = 0
   !====================================================================
   !lambda_minus = lambda0
   !====================================================================
   ! Procura um lambda_minus apropriado
   !do
   !   i1 = i1 + 1
   !   lambda_minus = lambda0
   !   !write(*,*) 'M Pot p achar lbd_-'
!=======================================================================
   !   call metodo_potencia(this, alp, mu, lambda_minus, autovalor_minus, tole)
!======================================================================= 
     !write(*,*) 'M pot concluiu' 
   !  if( autovalor_minus > 0.0d0 )then
   !      lambda0 = lambda0 - dlambda
   !      write(*,*) 'lbd_- = lbd_- - dlb'
   !   elseif( autovalor_minus == 0.0d0)then
   !      lambda_mean = lambda_minus
   !      return
   !   else
   !      write(*,*) "Achou lbd_-"
   !      exit
   !   endif
      
   !   if( lambda0 <= 0.0d0 )then
   !      dlambda = dlambda/2.0d0
   !      lambda0 = lambda_backup - dlambda
   !      write(*,*) 'Corrigindo lbd_- negativo'
   !   elseif( lambda0 >= 1.0d0)then
   !      lambda0 = lambda0 - dlambda0
   !      dlambda = dlambda0
   !   endif
   !enddo
!=======================================================================
   !i1 = 0
!=======================================================================

   !lambda0 = 2.0d0 * limiarQMF(this)
   
   !dlambda0 = lambda0/10.0d0
   !dlambda = dlambda0
     
   !lambda_backup = lambda0
   
   !====================================================================   
   !i1 = 0
   !====================================================================
   !lambda_plus = lambda0
   !====================================================================
   ! Procura um lambda_plus apropriado
   !do
   !   i1 = i1 + 1
      !write(*,*) 'M Pot para achar lbd_+'
!======================================================================= 
   !   call metodo_potencia(this, alp, mu, lambda_plus, autovalor_plus, tole)
!=======================================================================
   !   write(*,*) 'Concluiu M Pot'
   !   if( autovalor_plus < 0.0d0 )then
   !      write(*,*) 'lbd_+ = lbd_+ + dlb'
   !      lambda0 = lambda0 + dlambda
   !      dlambda = 2.0d0 * dlambda
   !   elseif( autovalor_plus == 0.0d0)then
   !      lambda_mean  = lambda_plus
   !      return
   !   else
   !      write(*,*) "Achou lbd_+"
   !      exit
   !   endif
    
   !   lambda_plus = lambda0
            
   !enddo
!=======================================================================
   !lambda_mean = (lambda_minus + lambda_plus)/2.0_dp
!=======================================================================
   !i1 = 0   
   
   lambda_mean = lambda0
   
   lambda_minus = 0.0d0
   
   lambda_plus = 0.0d0   
   !====================================================================
   l_biss_PQMF: do while(.True.)
   
      i1 = i1 + 1
      write(*,*) "Entrou no metodo da potencia"
!=======================================================================
     !call metodo_potencia(this, alp1, mu1, lambda1, autovalor1, tole)
      call metodo_potencia(this, alp, mu, lambda_mean, autovalor, tole)
!=======================================================================
      !write(*,*) "Saiu do metodo da potencia com autovalor = ", autovalor
      write(789, *) i1, lambda_mean

      if( ((lambda_minus * lambda_plus) == 0.0d0) )then
         if( autovalor < 0.0d0)then
            lambda_minus = lambda_mean
            autovalor_minus = autovalor
            lambda_mean = lambda_mean + dlambda
         else
            if( autovalor < tole )then
               write(*,*) 'Lambda_mean convergiu para = ', lambda_mean
               close(789)
               exit l_biss_PQMF
            else
               lambda_plus = lambda_mean
               autovalor_plus = autovalor
               dlambda = dlambda/2.0d0
               lambda_mean = lambda_mean - dlambda            
            endif
         endif                  
      else      
        if( (autovalor_minus * autovalor) < 0.0d0 )then
            lambda_plus = lambda_mean
            autovalor_plus = autovalor
         elseif( autovalor == 0.0d0 )then
            exit l_biss_PQMF
         else
            lambda_minus = lambda_mean
            autovalor_minus = autovalor
         endif
         lambda_mean = (lambda_minus + lambda_plus)/2.0d0
!=======================================================================
         if( abs( lambda_plus - lambda_minus ) < tole)then
            write(*,*) 'Lambda_mean convergiu para = ', lambda_mean
            close(789)
            exit l_biss_PQMF
         endif
!=======================================================================      
      endif
!=======================================================================      
   enddo l_biss_PQMF
   !====================================================================         
   if(allocated(fi)) deallocate(fi)
   allocate(fi(size(x)))
   !====================================================================         
   x_norm = sum(x**2.0d0)**0.5d0
   y_norm = sum(y**2.0d0)**0.5d0
   !====================================================================         
   fi = (y/y_norm + x/x_norm)/2.0d0
   deallocate(x)
   deallocate(y)
   write(*,*) "Convergiu com erro igual a ", abs(lambda_plus - lambda_minus)
   !====================================================================

   end subroutine

   function limiarQMF(this)
      !=================================================================
      class(grafo), intent(in) :: this
      real(dp) :: limiarQMF
      real(dp) :: P_grau(this%degMin:this%degMax)
      integer :: k1, k2, k3
      real(dp) :: qm, q2m
      
      !=================================================================
            
      !=================================================================      
      P_grau = 0.0_dp
      !=================================================================      
      do k1 = 1, this%nodes
         P_grau(this%deg(k1)) = P_grau(this%deg(k1)) + 1.0d0
      enddo
      !=================================================================      
      P_grau = P_grau/(1.0d0 * this%nodes)
      !=================================================================
      select type(this)
         class is (grafo_PL_UCM)         
            if( this%gamm > 2.5d0)then
               limiarQMF = 1.0d0/(this%degMax)**0.5d0
               write(*,*) "Usando como lambda0 limiar QMF = 1/sqrt(degMax) = ", limiarQMF
               
            else
               qm = 0.0d0
               q2m = 0.0d0
               do k1 = this%degMin, this%degMax
                  if(P_grau(k1) == 0.0d0) cycle
                  qm = qm + (1.0d0 * k1) * P_grau(k1)
                  q2m = q2m +  P_grau(k1) * (1.0d0 * k1)**2.0d0
               enddo
               limiarQMF = qm/q2m
               write(*,*) "Usando como lambda0 Limiar QMF = qm/q2m = ", limiarQMF
                                    
            endif
         class is (grafo_PL_UCM_Tubes)
            write(*,*) "Deu grafo PL UCM com Tubos"
            stop
         class is (grafo)
            write(*,*) "Deu classe grafo"
            stop
      end select
      !=================================================================
      
   end function
   !====================================================================
   subroutine aloca_listas_e_matrizes(this)
      class(grafo), intent(in) :: this
 
      if(allocated(x)) deallocate(x)
      
      allocate(x(this%nodes))

      if(allocated(y)) deallocate(y)
      
      allocate(y(this%nodes))      
 
   end subroutine
   
   !====================================================================
   subroutine metodo_potencia(this, alp1, mu1, lambda1, autovalor1, tole)
      !=================================================================
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp1, mu1, lambda1
      real(dp), intent(inout) :: autovalor1
      real(dp) :: tole
      !=================================================================
      real(dp) :: autovalor
      real(dp) :: autovalor_0      
      !=================================================================
      integer :: j1, j2, j3, j4
      integer(kind=8) :: j12
      real(dp) :: xp, yp
      integer :: ipx, ipy
      type(rndgen) :: geni
      real(dp) :: soma
      integer :: vizim
      real(dp) :: a_l, a_2l
      integer :: semente2
      real(dp) :: erro
      real(dp), allocatable :: v1(:)
      real(dp) :: x_norm, y_norm
      real(dp), allocatable :: x_prov(:)
      integer :: num_cont
      real(dp) :: Y4_0
      !=================================================================   
      a_l = func_a_l(alp1, lambda1, mu1)
      a_2l = func_a_2l(alp1, lambda1, mu1)      
      !=================================================================
      
      !=================================================================      
      call geni%init(semente2)     
      !=================================================================
      do j1 = 1, this%nodes
         y(j1) = this%deg(j1)
      enddo
      !=================================================================
      y_norm = (sum(y**2.0d0))**0.5d0
      !=================================================================
      if(allocated(x_prov)) deallocate(x_prov)
      allocate(x_prov(this%nodes))
      !=================================================================
      ! Apos normalizar o vetor x, a norma ||xp|| = ||x||_{oo}
      ! se torna = 1.0d0.
      !=================================================================
      j1 = 0
      !=================================================================
      num_cont = 2
      Y4_0 = 0.0d0
      autovalor_0 = 0.0d0
!=======================================================================      
iter: do while(.True.)	
	     !Aqui comeca o algoritmo
         !==============================================================
         ! Y = A X
         !==============================================================
         y = y/y_norm
         x_prov = y
         !==============================================================
         do j3 = 1, num_cont
            x = y
            do j1 = 1, this%nodes
               !========================================================
               soma = 0.0_dp
               do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
                  vizim = this%listAdj(j12)
                  soma = soma + x(vizim)
               enddo
               y(j1) = a_l * soma
               !========================================================
               y(j1) = y(j1) + a_2l * ( 1.0d0 * ( this%degMax - this%deg(j1)) ) * x(j1)
               !========================================================
            enddo
            
         enddo
         !==============================================================         
         y_norm = (sum(y**2.0d0))**0.5d0
         x_norm = (sum(x**2.0d0))**0.5d0
         Y4 = sum(((y/y_norm + x/x_norm)/2.0d0)**4.0d0)
         !==============================================================
         autovalor = (dot_product(x_prov, y))**(1.0d0/(1.0d0 * num_cont))
         !==============================================================
         erro = max(abs(autovalor - autovalor_0), abs(Y4 - Y4_0))
         autovalor_0 = autovalor
         Y4_0 = Y4
         !==============================================================                  
         if(erro < tole)then
            autovalor1 = autovalor - ( a_2l * (1.0d0 * this%degMax) + mu1 )
            write(*,*) "Eig1 = ", autovalor1         
            !-----------------------------------------------------------
            exit iter
         endif         
	enddo iter      
      
      
   end subroutine

   function func_a_l(alfa, lambda, mu)
      real(dp), intent(in) :: alfa
      real(dp), intent(in) :: lambda
      real(dp), intent(in) :: mu
      real(dp) :: func_a_l
      real(dp) :: a_par
      
      !=======================================================================
      a_par = (alfa * mu)/(alfa + 2.0_dp * mu)
      !=======================================================================
             
      func_a_l = lambda * ( ( a_par + lambda + mu )/( a_par + 2.0_dp * lambda + mu ) )      
          
   end function
   
   function func_a_2l(alfa, lambda, mu)
      real(dp), intent(in) :: alfa
      real(dp), intent(in) :: lambda
      real(dp), intent(in) :: mu
      real(dp) :: func_a_2l
      real(dp) :: a_par
      
      !=======================================================================
      a_par = (alfa * mu)/(alfa + 2.0_dp * mu)
      !=======================================================================
                   
      func_a_2l = (lambda ** 2.0_dp)/( a_par + 2.0_dp * lambda + mu )      
      
   end function   
 
end module

module sirs_pqmf
   use geraRede
   use mod_rndgen
   use mod_tools
   implicit none
   
   integer, parameter :: dp = kind(0.0d0)
   
   real(dp), allocatable :: Ii(:), Ri(:)
   real(dp), allocatable :: Ii_tild(:), Ri_tild(:)
   real(dp), allocatable :: SIij(:), RRij(:), RIij(:)
   real(dp), allocatable :: RIij_tild(:), RRij_tild(:), SIij_tild(:)
   integer, allocatable :: spec(:)
   
!#######################################################################  
   real(dp), allocatable :: k1_Ii(:), k2_Ii(:), k3_Ii(:), k4_Ii(:)
   real(dp), allocatable :: k1_Ri(:), k2_Ri(:), k3_Ri(:), k4_Ri(:)
!#######################################################################      
   real(dp), allocatable :: k1_RIij(:), k2_RIij(:), k3_RIij(:), k4_RIij(:)
   real(dp), allocatable :: k1_SIij(:), k2_SIij(:), k3_SIij(:), k4_SIij(:)
   real(dp), allocatable :: k1_RRij(:), k2_RRij(:), k3_RRij(:), k4_RRij(:)
!#######################################################################   
   real(dp), allocatable :: v1(:)
   
   real(dp) :: t, dt
   
   contains

!#######################################################################   
   subroutine aloca_listas(this)
      class(grafo), intent(in) :: this
      
      integer :: n_sitios
      integer (kind=8) :: n_arestas
      integer :: k1, k2, k3
      integer(kind=8) :: k12, k13, k23
      
      n_sitios = this%nodes
      
      n_arestas = this%sumDeg
      
      !#################################################################
      if(allocated(Ii)) deallocate(Ii)
         allocate(Ii(n_sitios))

      if(allocated(Ri)) deallocate(Ri)
         allocate(Ri(n_sitios))         

      if(allocated(SIij)) deallocate(SIij)
         allocate(SIij(n_arestas))

      if(allocated(RRij)) deallocate(RRij)
         allocate(RRij(n_arestas))         

      if(allocated(RIij)) deallocate(RIij)
         allocate(RIij(n_arestas))   
      
      if(allocated(spec)) deallocate(spec)
         allocate(spec(n_arestas))
      
      do k1 = 1, this%nodes
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            k2 = this%listAdj(k12)
l_k23:      do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               k3 = this%listAdj(k23)
               if( k3 == k1 )then
                  spec(k12) = k23
                  exit l_k23                  
               endif
            enddo l_k23
         enddo 
      enddo

      if(allocated(k1_Ii))then
         deallocate(k1_Ii) 
         deallocate(k2_Ii)
         deallocate(k3_Ii)
         deallocate(k4_Ii)
      endif
      allocate(k1_Ii(this%nodes))
      allocate(k2_Ii(this%nodes))
      allocate(k3_Ii(this%nodes))
      allocate(k4_Ii(this%nodes))
      

      if(allocated(k1_Ri))then
         deallocate(k1_Ri) 
         deallocate(k2_Ri)
         deallocate(k3_Ri)
         deallocate(k4_Ri)
      endif
      allocate(k1_Ri(this%nodes))
      allocate(k2_Ri(this%nodes))
      allocate(k3_Ri(this%nodes))
      allocate(k4_Ri(this%nodes))
      

      if(allocated(k1_RIij))then
         deallocate(k1_RIij)
         deallocate(k2_RIij)
         deallocate(k3_RIij)
         deallocate(k4_RIij)
      endif
      allocate(k1_RIij(n_arestas))
      allocate(k2_RIij(n_arestas))
      allocate(k3_RIij(n_arestas))
      allocate(k4_RIij(n_arestas))

      if(allocated(k1_RRij))then
         deallocate(k1_RRij)
         deallocate(k2_RRij)
         deallocate(k3_RRij)
         deallocate(k4_RRij)
      endif 
      allocate(k1_RRij(n_arestas))
      allocate(k2_RRij(n_arestas))
      allocate(k3_RRij(n_arestas))
      allocate(k4_RRij(n_arestas))      

      if(allocated(k1_SIij))then
         deallocate(k1_SIij)
         deallocate(k2_SIij)
         deallocate(k3_SIij)
         deallocate(k4_SIij)
      endif
      allocate(k1_SIij(n_arestas))
      allocate(k2_SIij(n_arestas))
      allocate(k3_SIij(n_arestas))
      allocate(k4_SIij(n_arestas))
      !=================================================================            
   end subroutine
   
   subroutine condicao_inicial(this, I_0, t1)     

      class(grafo), intent(in) :: this
      real(dp), intent(in) :: I_0
      integer:: k1, k2
      integer(kind=8) :: k12, k21, k23
      logical :: isI_i
      integer :: iost
      character(len=30) :: Nichar
      character(len=10) :: lambdachar
      character(len=100) :: arquivo_char2
      integer(kind = 8) :: sum_deg
      real(dp) :: t1    
      t = t1
                     
      do k1 = 1, this%nodes
         if(lista_de_clusters(k1) /= i_comp_gigante)then
            Ii(k1) = 0.0_dp
            Ri(k1) = 0.0_dp
            cycle
         endif
         Ii(k1) = I_0
         Ri(k1) = 0.01_dp
      enddo

      do k1 = 1, this%nodes

         if(lista_de_clusters(k1) /= i_comp_gigante) cycle

         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            
            k2 = this%listAdj(k12)
                                         
            SIij(k12) = (1.0_dp - Ii(k1) - Ri(k1)) * Ii(k2)            
               
            RIij(k12) = Ri(k1) * Ii(k2)!0.0_dp
 
            RRij(k12) = Ri(k1) * Ri(k2) !0.0_dp
            
         enddo            

      enddo
                               
   end subroutine    
   
   
!#######################################################################   
   subroutine rk4(dt, t, this, n_sitios, n_arestas, alf, lbd, mu)
      real(dp) :: dt
      real(dp) :: t
      class(grafo), intent(in) :: this         
      integer :: n_sitios
      integer (kind=8) :: n_arestas
      real(dp) :: alf, lbd, mu
      real(dp), parameter :: factor = 1.0d0/6.0d0
      
         !#########################SUBROTINA#############################

!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, SIij)

		k1_Ii = dt * f_Ii(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, SIij)
		
		
		            !f_Ri(this, n_sitios, t, alf, mu, Ii, Ri)

		k1_Ri = dt * f_Ri(this, n_sitios, t, alf, mu, Ii, Ri)

!#######################################################################		
                      !f_RIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)

		k1_RIij = dt * f_RIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)

!#######################################################################		
                      !f_RRij(this, n_sitios, n_arestas, t, alf, mu, RIij, RRij)

		k1_RRij = dt * f_RRij(this, n_sitios, n_arestas, t, alf, mu, RIij, RRij)

!#######################################################################		
                      !f_SIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)

		k1_SIij = dt * f_SIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)


!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t,               alf, lbd, mu, Ii                 , SIij                   )

		k2_Ii = dt * f_Ii(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k1_Ii, SIij + 0.5_dp * k1_SIij)

!#######################################################################
                    !f_Ri(this, n_sitios, t,               alf, mu, Ii                 , Ri                 )

		k2_Ri = dt * f_Ri(this, n_sitios, t + 0.5_dp * dt, alf, mu, Ii + 0.5_dp * k1_Ii, Ri + 0.5_dp * k1_Ri)

!#######################################################################		
                      !f_RIij(this, n_sitios, n_arestas, t,               alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )

		k2_RIij = dt * f_RIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k1_Ii, Ri + 0.5_dp * k1_Ri, RIij + 0.5_dp * k1_RIij, RRij + 0.5_dp * k1_RRij, SIij + 0.5_dp * k1_SIij)

!#######################################################################
                      !f_RRij(this, n_sitios, n_arestas, t,               alf, mu, RIij,                    RRij                   )

		k2_RRij = dt * f_RRij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, mu, RIij + 0.5_dp * k1_RIij, RRij + 0.5_dp * k1_RRij)


!#######################################################################		
                      !f_SIij(this, n_sitios, n_arestas, t              , alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )

		k2_SIij = dt * f_SIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k1_Ii, Ri + 0.5_dp * k1_Ri, RIij + 0.5_dp * k1_RIij, RRij + 0.5_dp * k1_RRij, SIij + 0.5_dp * k1_SIij)


!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t             , alf, lbd, mu,  Ii                 , SIij                   )
		
		k3_Ii = dt * f_Ii(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k2_Ii, SIij + 0.5_dp * k2_SIij)

!#######################################################################
		             !f_Ri(this, n_sitios, t             , alf, mu, Ii                 , Ri                 )
		
		k3_Ri = dt * f_Ri(this, n_sitios, t + 0.5_dp * dt, alf, mu, Ii + 0.5_dp * k2_Ii, Ri + 0.5_dp * k2_Ri)
		
!#######################################################################
                      !f_RIij(this, n_sitios, n_arestas, t              , alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )

		k3_RIij = dt * f_RIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k2_Ii, Ri + 0.5_dp * k2_Ri, RIij + 0.5_dp * k2_RIij, RRij + 0.5_dp * k2_RRij, SIij + 0.5_dp * k2_SIij)

!#######################################################################
                      !f_RRij(this, n_sitios, n_arestas, t              , alf, mu, RIij                   , RRij                   )

		k3_RRij = dt * f_RRij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, mu, RIij + 0.5_dp * k2_RIij, RRij + 0.5_dp * k2_RRij)

!#######################################################################
                      !f_SIij(this, n_sitios, n_arestas, t              , alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )
        
		k3_SIij = dt * f_SIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k2_Ii, Ri + 0.5_dp * k2_Ri, RIij + 0.5_dp * k2_RIij, RRij + 0.5_dp * k2_RRij, SIij + 0.5_dp * k2_SIij)

!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t     , alf, lbd, mu, Ii        , SIij          )
		k4_Ii = dt * f_Ii(this, n_sitios, n_arestas, t + dt, alf, lbd, mu, Ii + k3_Ii, SIij + k3_SIij)
!#######################################################################						
                    !f_Ri(this, n_sitios, t     , alf, mu, Ii        , Ri        )

		k4_Ri = dt * f_Ri(this, n_sitios, t + dt, alf, mu, Ii + k3_Ii, Ri + k3_Ri)

!#######################################################################		
                      !f_RIij(this, n_sitios, n_arestas, t     , alf, lbd, mu, Ii        , Ri        , RIij          , RRij          , SIij          )

		k4_RIij = dt * f_RIij(this, n_sitios, n_arestas, t + dt, alf, lbd, mu, Ii + k3_Ii, Ri + k3_Ri, RIij + k3_RIij, RRij + k3_RRij, SIij + k3_SIij)


!#######################################################################				                
                      !f_RRij(this, n_sitios, n_arestas, t     , alf, mu, RIij          , RRij          )

		k4_RRij = dt * f_RRij(this, n_sitios, n_arestas, t + dt, alf, mu, RIij + k3_RIij, RRij + k3_RRij)


!#######################################################################	
                      !f_SIij(this, n_sitios, n_arestas, t     , alf, lbd, mu, Ii        , Ri        , RIij          , RRij          , SIij          )
        
		k4_SIij = dt * f_SIij(this, n_sitios, n_arestas, t + dt, alf, lbd, mu, Ii + k3_Ii, Ri + k3_Ri, RIij + k3_RIij, RRij + k3_RRij, SIij + k3_SIij)
						
!################################################################################################
! Factor = 1.0d0/6.0d0 - it is a parameter
!################################################################################################
        Ii =   Ii   + factor * (k1_Ii   + 2.0_dp * k2_Ii   + 2.0_dp * k3_Ii   + k4_Ii)

        Ri =   Ri   + factor * (k1_Ri   + 2.0_dp * k2_Ri   + 2.0_dp * k3_Ri   + k4_Ri)

        RIij = RIij + factor * (k1_RIij + 2.0_dp * k2_RIij + 2.0_dp * k3_RIij + k4_RIij)
        
        RRij = RRij + factor * (k1_RRij + 2.0_dp * k2_RRij + 2.0_dp * k3_RRij + k4_RRij)
        
        SIij = SIij + factor * (k1_SIij + 2.0_dp * k2_SIij + 2.0_dp * k3_SIij + k4_SIij)         
!#########################SUBROTINA##############################################################      
   end subroutine
   !####################################################################
   !   Funcoes
   !####################################################################
   subroutine rk4_teste(dt, t, this, alf, lbd, mu)
      real(dp) :: dt
      real(dp) :: t
      class(grafo), intent(in) :: this         
      integer :: n_sitios
      integer (kind=8) :: n_arestas
      real(dp) :: alf, lbd, mu
      real(dp), parameter :: factor = 1.0d0/6.0d0
      !#################################################################
      n_sitios = this%nodes; n_arestas = this%sumDeg
!#######################################################################
                    !f_Ii(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, SIij)

		k1_Ii = dt * f_Ii(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, SIij)
		            !f_Ri(this, n_sitios, t, alf, mu, Ii, Ri)

		k1_Ri = dt * f_Ri(this, n_sitios, t, alf, mu, Ii, Ri)

!#######################################################################		
                      !f_RIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)

		k1_RIij = dt * f_RIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)

!#######################################################################		
                      !f_RRij(this, n_sitios, n_arestas, t, alf, mu, RIij, RRij)

		k1_RRij = dt * f_RRij(this, n_sitios, n_arestas, t, alf, mu, RIij, RRij)

!#######################################################################		
                      !f_SIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)

		k1_SIij = dt * f_SIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)


!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t,               alf, lbd, mu, Ii                 , SIij                   )

		k2_Ii = dt * f_Ii(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k1_Ii, SIij + 0.5_dp * k1_SIij)

!#######################################################################
                    !f_Ri(this, n_sitios, t,               alf, mu, Ii                 , Ri                 )

		k2_Ri = dt * f_Ri(this, n_sitios, t + 0.5_dp * dt, alf, mu, Ii + 0.5_dp * k1_Ii, Ri + 0.5_dp * k1_Ri)

!#######################################################################		
                      !f_RIij(this, n_sitios, n_arestas, t,               alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )

		k2_RIij = dt * f_RIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k1_Ii, Ri + 0.5_dp * k1_Ri, RIij + 0.5_dp * k1_RIij, RRij + 0.5_dp * k1_RRij, SIij + 0.5_dp * k1_SIij)

!#######################################################################
                      !f_RRij(this, n_sitios, n_arestas, t,               alf, mu, RIij,                    RRij                   )

		k2_RRij = dt * f_RRij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, mu, RIij + 0.5_dp * k1_RIij, RRij + 0.5_dp * k1_RRij)


!#######################################################################		
                      !f_SIij(this, n_sitios, n_arestas, t              , alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )

		k2_SIij = dt * f_SIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k1_Ii, Ri + 0.5_dp * k1_Ri, RIij + 0.5_dp * k1_RIij, RRij + 0.5_dp * k1_RRij, SIij + 0.5_dp * k1_SIij)


!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t             , alf, lbd, mu,  Ii                 , SIij                   )
		
		k3_Ii = dt * f_Ii(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k2_Ii, SIij + 0.5_dp * k2_SIij)

!#######################################################################
		             !f_Ri(this, n_sitios, t             , alf, mu, Ii                 , Ri                 )
		
		k3_Ri = dt * f_Ri(this, n_sitios, t + 0.5_dp * dt, alf, mu, Ii + 0.5_dp * k2_Ii, Ri + 0.5_dp * k2_Ri)
		
!#######################################################################
                      !f_RIij(this, n_sitios, n_arestas, t              , alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )

		k3_RIij = dt * f_RIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k2_Ii, Ri + 0.5_dp * k2_Ri, RIij + 0.5_dp * k2_RIij, RRij + 0.5_dp * k2_RRij, SIij + 0.5_dp * k2_SIij)

!#######################################################################
                      !f_RRij(this, n_sitios, n_arestas, t              , alf, mu, RIij                   , RRij                   )

		k3_RRij = dt * f_RRij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, mu, RIij + 0.5_dp * k2_RIij, RRij + 0.5_dp * k2_RRij)

!#######################################################################
                      !f_SIij(this, n_sitios, n_arestas, t              , alf, lbd, mu, Ii                 , Ri                 , RIij                   , RRij                   , SIij                   )
        
		k3_SIij = dt * f_SIij(this, n_sitios, n_arestas, t + 0.5_dp * dt, alf, lbd, mu, Ii + 0.5_dp * k2_Ii, Ri + 0.5_dp * k2_Ri, RIij + 0.5_dp * k2_RIij, RRij + 0.5_dp * k2_RRij, SIij + 0.5_dp * k2_SIij)

!#################################################################################################
                    !f_Ii(this, n_sitios, n_arestas, t     , alf, lbd, mu, Ii        , SIij          )
		k4_Ii = dt * f_Ii(this, n_sitios, n_arestas, t + dt, alf, lbd, mu, Ii + k3_Ii, SIij + k3_SIij)
!#######################################################################						
                    !f_Ri(this, n_sitios, t     , alf, mu, Ii        , Ri        )

		k4_Ri = dt * f_Ri(this, n_sitios, t + dt, alf, mu, Ii + k3_Ii, Ri + k3_Ri)

!#######################################################################		
                      !f_RIij(this, n_sitios, n_arestas, t     , alf, lbd, mu, Ii        , Ri        , RIij          , RRij          , SIij          )

		k4_RIij = dt * f_RIij(this, n_sitios, n_arestas, t + dt, alf, lbd, mu, Ii + k3_Ii, Ri + k3_Ri, RIij + k3_RIij, RRij + k3_RRij, SIij + k3_SIij)


!#######################################################################				                
                      !f_RRij(this, n_sitios, n_arestas, t     , alf, mu, RIij          , RRij          )

		k4_RRij = dt * f_RRij(this, n_sitios, n_arestas, t + dt, alf, mu, RIij + k3_RIij, RRij + k3_RRij)


!#######################################################################	
                      !f_SIij(this, n_sitios, n_arestas, t     , alf, lbd, mu, Ii        , Ri        , RIij          , RRij          , SIij          )
        
		k4_SIij = dt * f_SIij(this, n_sitios, n_arestas, t + dt, alf, lbd, mu, Ii + k3_Ii, Ri + k3_Ri, RIij + k3_RIij, RRij + k3_RRij, SIij + k3_SIij)
!################################################################################################
! Factor = 1.0d0/6.0d0 - it is a parameter
!################################################################################################
        Ii_tild =   Ii   + factor * (k1_Ii   + 2.0_dp * k2_Ii   + 2.0_dp * k3_Ii   + k4_Ii)

        Ri_tild =   Ri   + factor * (k1_Ri   + 2.0_dp * k2_Ri   + 2.0_dp * k3_Ri   + k4_Ri)

        RIij_tild = RIij + factor * (k1_RIij + 2.0_dp * k2_RIij + 2.0_dp * k3_RIij + k4_RIij)
        
        RRij_tild = RRij + factor * (k1_RRij + 2.0_dp * k2_RRij + 2.0_dp * k3_RRij + k4_RRij)
        
        SIij_tild = SIij + factor * (k1_SIij + 2.0_dp * k2_SIij + 2.0_dp * k3_SIij + k4_SIij)         
!#########################SUBROTINA##############################################################      
   end subroutine
   !####################################################################
   !   Funcoes
!#######################################################################   
   function f_Ii(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, SIij)
      class(grafo), intent(in) :: this      
      integer, intent(in) :: n_sitios
      integer (kind=8), intent(in) :: n_arestas
      real(dp) :: f_Ii(n_sitios)
      real(dp), intent(in) :: t
      real(dp), intent(in) :: alf
      real(dp), intent(in) :: lbd
      real(dp), intent(in) :: mu
      real(dp) :: Ii(n_sitios)
      real(dp) :: SIij(n_arestas)      
      integer :: k1
      integer(kind = 8) :: k12
      real(dp) :: SIij_aux
      
      !#################################################################            
      do k1 = 1, this%nodes
         !##############################################################
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
         !##############################################################
         SIij_aux = 0.0_dp
         !##############################################################
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            SIij_aux = SIij_aux + SIij(k12)
         enddo
         !##############################################################
         f_Ii(k1) = -mu * Ii(k1) + lbd * SIij_aux
         !##############################################################
      enddo
      !#################################################################      
   end function

!#######################################################################
   function f_Ri(this, n_sitios, t, alf, mu, Ii, Ri)
      class(grafo), intent(in) :: this
      integer, intent(in) :: n_sitios
      real(dp) :: f_Ri(n_sitios)
      real(dp), intent(in) :: t
      real(dp), intent(in) :: alf
      real(dp), intent(in) :: mu
      real(dp) :: Ii(n_sitios)
      real(dp) :: Ri(n_sitios)
      integer :: k1
      !#################################################################      
      do k1 = 1, this%nodes
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
         f_Ri(k1) = mu * Ii(k1) - alf * Ri(k1)
      enddo
      !#################################################################
   end function

!#######################################################################
   function f_SIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)
      class(grafo), intent(in) :: this      
      integer, intent(in) :: n_sitios
      integer (kind=8), intent(in) :: n_arestas
      real(dp) :: f_SIij(n_arestas)
      real(dp), intent(in) :: t
      real(dp), intent(in) :: alf
      real(dp), intent(in) :: lbd
      real(dp), intent(in) :: mu
      real(dp) :: Ii(n_sitios)
      real(dp) :: Ri(n_sitios)
      real(dp) :: RIij(n_arestas)
      real(dp) :: RRij(n_arestas)
      real(dp) :: SIij(n_arestas)      
      integer :: k1, k2, k3
      integer(kind=8) :: k12, k21, k13, k31, k23, k32
      real(dp) :: inv_S1, inv_S2
      
      real(dp) :: SIik_aux, SIjk_aux
      real(dp) :: SSij_aux
      
      do k1 = 1, this%nodes
         !##############################################################
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
         !##############################################################
         !##############################################################         
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            !###########################################################
            f_SIij(k12) = alf * RIij(k12) - ( lbd + mu ) * SIij(k12)
            !###########################################################
 
            k21 = spec(k12)
            
            k2 = this%listAdj(k12)
            !###########################################################
            SSij_aux = (1.0_dp - Ii(k1) - Ri(k1)) - SIij(k12) - ( Ri(k2) - RIij(k21) - RRij(k21) )
            ! SSij   = Si                         - SIij      -SRij
            
            ! Rj   = RIji + RRji + RSji
            ! RSji = Rj   - RIji - RRji
            ! SRij = Rj   - RIji - RRij
            
            ! RIij ===> RIij(k12)
            ! RIji ===> RIij(k21)
            ! ISij ===> SIji ===> SI(ji)
            !###########################################################
            SIjk_aux = 0.0_dp
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               if( k23 == k21 ) cycle
               SIjk_aux = SIjk_aux + SIij(k23)
            enddo
            !###########################################################
            f_SIij(k12) = f_SIij(k12) + lbd * SSij_aux * SIjk_aux/ ( 1.0_dp - Ii(k2) - Ri(k2) ) 
            !###########################################################
            SIik_aux = 0.0_dp
            do k13 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
               if( k13 == k12 ) cycle
               SIik_aux = SIik_aux + SIij(k13)
            enddo
            !###########################################################
            ! Antes tinha um sinal de -. Optei por mudar o sinal de Si.
            f_SIij(k12) = f_SIij(k12) - lbd * SIik_aux * SIij(k12) /( 1.0_dp - Ri(k1) - Ii(k1) ) 
            !###########################################################
         enddo
      enddo
   end function
!#######################################################################   
   function f_RRij(this, n_sitios, n_arestas, t, alf, mu, RIij, RRij)
      class(grafo), intent(in) :: this      
      integer, intent(in) :: n_sitios
      integer (kind=8), intent(in) :: n_arestas
      real(dp) :: f_RRij(n_arestas)
      real(dp), intent(in) :: t
      real(dp), intent(in) :: alf
      real(dp), intent(in) :: mu
      real(dp) :: RIij(n_arestas)
      real(dp) :: RRij(n_arestas)      
      integer :: k1, k2
      integer(kind=8) :: k12, k21
      
      !#################################################################
      do k1 = 1, this%nodes
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            k21 = spec(k12)
            f_RRij(k12) = mu * (RIij(k21) + RIij(k12)) -2.0_dp  * alf * RRij(k12)
         enddo
      enddo
      !#################################################################   
   end function
!#######################################################################
   function f_RIij(this, n_sitios, n_arestas, t, alf, lbd, mu, Ii, Ri, RIij, RRij, SIij)
      class(grafo), intent(in) :: this      
      integer, intent(in) :: n_sitios
      integer (kind=8), intent(in) :: n_arestas
      real(dp) :: f_RIij(n_arestas)
      real(dp), intent(in) :: t
      real(dp), intent(in) :: alf
      real(dp), intent(in) :: lbd
      real(dp), intent(in) :: mu
      real(dp) :: Ii(n_sitios)
      real(dp) :: Ri(n_sitios)
      real(dp) :: RIij(n_arestas)
      real(dp) :: RRij(n_arestas)
      real(dp) :: SIij(n_arestas)      
      integer :: k1, k2, k3
      integer(kind = 8) :: k12, k21, k13, k31, k23, k32
      
      real(dp) :: SIjk_aux
      real(dp) :: SSij_aux
      
      !open(444, file=trim(adjustl(local))//'absIIij_IIji.dat', status='unknown')
      
      do k1 = 1, this%nodes
         !###########################################################
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
         !###########################################################
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            !###########################################################            
            k2 = this%listAdj(k12)
            k21 = spec(k12)
            !###########################################################
            f_RIij(k12) = mu * ( Ii(k1) - RIij(k21) - SIij(k21) ) -( alf + mu ) * RIij(k12) 
            !###########################################################
            SIjk_aux = 0.0_dp
            !###########################################################
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               if( k23 == k21 ) cycle               
               SIjk_aux = SIjk_aux + SIij(k23)               
            enddo
            !###########################################################
            f_RIij(k12) = f_RIij(k12) + lbd  * SIjk_aux * ( Ri(k1) -RIij(k12) - RRij(k12) )/( 1.0_dp - Ii(k2) - Ri(k2) )
            !###########################################################
         enddo
      enddo
      !close(444)
!#######################################################################
   end function
end module

program main
   use geraRede
   use sirs_pqmf
   use mod_rndgen
   use mod_tools
   use limiaresCampoMedio
   implicit none
!#######################################################################   
   type(grafo_PL_UCM) :: rede
!#######################################################################   
   real(dp) :: dlamb
   real(dp) :: lamb0, lambdaf
   real(dp) :: lamb
   integer :: nlamb = 1000
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
   real(dp) :: dt0
   real(dp) :: rho, rho0, reco
   real(dp) :: tol, tole
!#######################################################################
   !type(rndgen) :: gen1
   integer :: seed(50)
   integer :: seed1
   type(rndgen) :: ger_inic
!#######################################################################
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   real(dp) :: grau_max      
!#######################################################################   
   character(len=500) :: t_vs_Im
   character(len=500) :: lamb_vs_Im
   character(len=500) :: lamb_vs_Xi 
   character(len=1000) :: caminho
   character(len=500) :: arquivo_rede_real
   character(len=5) :: tam_char
   character(len=5) :: gama_char
   character(len=5) :: indice
!#######################################################################   
   integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9 !, j1, j2, j3
   integer(kind = 8) :: i12, i21, i13, i31, i23, i32
   logical :: T_vs
!#######################################################################
   integer :: sumDeg2
   real(dp) :: t0, tf
   real(dp) :: dt_m
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   integer :: ntempo
   
   integer :: n_it_rho, n_it_Y4, n_sub, per_conv
   
   real(dp) :: p0
   character(len=300) :: cwd, resultados, tipoCorte
   character(len=1000) :: local
   character(len=1000) :: nomeArquivo
   character(len=20) :: buffer
   !##############################################################
   
   real(dp) :: qm, q2m

   integer :: nargus, ind_lamb   
   character(len=10) :: char_ind_lamb
   integer :: niter_abs   

   real(dp) :: lbdC_PQMF
   integer :: ind_amostra
   integer :: divisor
   character(len=10) :: lamb_char
   real(dp) ::lamb_C
   real(dp) :: tempoEscrita, tempoEscrita0
   real(dp) :: lamb_Copia
   logical :: existe
   integer :: st
   logical :: usouCopia
   logical :: teveLeitura
   real(dp) :: lambda_Ultimo_Index
   character(len=10) :: teoriaCM

   integer :: int_soCalculaIPR
   logical :: soCalculaIPR
   logical :: taAberta
   real(dp) :: Y4_old
   real(dp) :: sumIi2
   logical :: Ii_neg
   logical :: Ii_gt_1
   integer(kind=8) :: so_gasta
   real(dp) :: lbd_lido, Y4_lido
   character(len = 100) :: C_lamb
   character(len = 100) :: C_tempo, C_tempoEscrita
   character(len = 100) :: C_I_i, C_R_i, C_RRij, C_RIij, C_SIij
   integer :: st10
   logical :: existe10
   integer :: gasta_aleatorio
   real(dp) :: vec_dt(50)
   integer :: ind_degMax
   real(dp) :: Ii_tilda
   real(dp) :: delta_Ii
   real(dp) :: claus_Ii
!=======================================================================
! ----------------Name space dos arquivos de backup---------------------
!=======================================================================
   C_lamb = 'Copia_lambda'
   C_tempo = 'Copia_tempo'
   C_I_i = 'Copia_I_i'
   C_R_i = 'Copia_R_i'
   C_RRij = 'Copia_RR_ij'
   C_RIij = 'Copia_RI_ij'
   C_SIij = 'Copia_SI_ij'
   C_tempoEscrita = 'Copia_tempoEscrita'
!=======================================================================
   seed1=947361823
!=======================================================================
   salvouCentralidade = .False.
!=======================================================================
 
   teoriaCM = trim(adjustl('pQMF'))
   
   !tipoCorte ='_Rigido'
   tipoCorte ='_2sqrtN'
   !tipoCorte ='_sqrtN'
   
   resultados = 'Rst_'//trim(adjustl(teoriaCM))//'_Corte'//trim(adjustl(tipoCorte))
    
   resultados = trim(adjustl(resultados))

   call system('mkdir -p '//resultados)
    
   local = trim(adjustl(resultados))//"/"


   call entradaArgumentos()

!#######################################################################
   !====================================================================
   !  Se der ruim no dt, ele eh dividido por dois.
   !====================================================================   
   dt0 = 1.0d-1 
   vec_dt(1) = dt0
   do i1 = 2, size(vec_dt)
      vec_dt(i1) = vec_dt(i1 - 1) * 0.5d0
   enddo
   !====================================================================
   ! A principio, t0 = 0, mas, se houver um estado salvo, muda
   ! para o t salvo no arquivo.
   !====================================================================
   t0 = 0.0_dp
   tf = 10000.0_dp
   
   tempoEscrita0 = 5.0d0
   
   ntempo = 100 * int( (tf- t0)/dt0 )
   tole = 1d-7
   per_conv = int( 10.0_dp/dt0 )
!#######################################################################

!#######################################################################
!   call ger_inic%init(seed1)
!   i2 = 1
!   do i1 = 1, 1000
!      if(mod(i1,100) > 0) cycle
!      seed(i2)  = ger_inic%int(100000000,999999999)
!      write(*,*) i1, seed(i2)
!      i2 = i2+1      
!   enddo
!   i2 = 1
!   do i1 = 1, 1000
!      if(mod(i1,100) /= 50)then !Antes, eu dava um cycle. Agora eh assim      
!         so_gasta = ger_inic%int(100000000,999999999)
!      else
!         seed(i2)  = ger_inic%int(100000000,999999999)
!         write(*,*) i1, seed(i2)
!         i2 = i2+1
!      endif
!   enddo
   !====================================================================
   call ger_inic%init(seed1)
   !====================================================================   
   if(.True.)then
      !=================================================================
      ! Vamos usar esse,
      ! que parece ser a versao antiga (mas que pode concordar com algumas
      ! amostras de rede que temos.
      !=================================================================
      i2 = 1
      do i1 = 1, 5000
         if(mod(i1,100) > 0) cycle
         seed(i2)  = ger_inic%int(100000000,999999999)
         write(*,*) i1, seed(i2)
         i2 = i2+1      
      enddo
   else
      !=================================================================
      ! Ou esse,
      ! que eh da versao mais nova e que concorda com algumas redes.
      !=================================================================   
      i2 = 1
      do i1 = 1, 1000
         if(mod(i1,100) > 0)then
            gasta_aleatorio = ger_inic%int(100000000,999999999)
            cycle
         endif
         seed(i2)  = ger_inic%int(100000000,999999999)
         write(*,*) i1, seed(i2)
         i2 = i2+1
      enddo
   endif   
   !====================================================================
   local = trim(adjustl(trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
   
   local = trim(adjustl(trim(adjustl(local))//'gam_'//trim(adjustl(gama_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )   

   local = trim(adjustl(trim(adjustl(local))//'ams_'//trim(adjustl(indice))//'/'))
   
   call system('mkdir -p '//trim(adjustl(local)) )
!#######################################################################
!				Inicia grafo
!#######################################################################
 !  arquivo_rede_real='s00088.s838.net.edg'
 !  call rede%RedeReal(arquivo_rede_real, 111)
 !  close(111)
!#######################################################################
   call criaRedeEClassificaClusters(rede, tam_rede, grau_min, grau_max, gama_exp, seed(ind_amostra))

   acha_hub:do i1 = 1, rede%nodes
      if(rede%deg(i1) == rede%degMax)then
         ind_degMax = i1
         exit acha_hub
      endif
   enddo acha_hub
   !====================================================================
   write(*,*) "######################Dados da Rede######################"
   write(*,*) ""
   write(*,*) "Tamanho da rede ", rede%nodes, "."
   write(*,*) ""
   write(*,*) "Fracao correspondente aa componente gigante ", 100.0 * comp_gigante/rede%nodes,"%", "."
   write(*,*) ""
   write(*,*) "Expoente da distribuicao de graus da rede ", gama_exp, "."
   write(*,*) ""
   write(*,*) "Grau minimo ", rede%degMin, ".", " Grau maximo ", rede%degMax, "."
   write(*,*) ""
!#######################################################################
   write(alp_char2, '(f9.3)') alp

   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//trim(adjustl(tipoCorte))//'/'))
 
   call system('mkdir -p '//trim(adjustl(local)) )
  
   call kNN_e_clustering(rede)

!#######################################################################
   write(*,*) "######################Dados temporais######################"
   write(*,*) ""
   write(*,*) "Instante inicial ", t0, "."	
   write(*,*) ""
   write(*,*) "Resolucao ", dt, "."
   write(*,*) ""
   write(*,*) "Instante final ", tf, "."
   write(*,*) ""
   write(*,*) "Quantidade de pontos ", ntempo, "."
   write(*,*) ""
   write(*,*) "Tolerancia considerada no criterio de conv. ", tole, "."
!#######################################################################
   write(*,*) "###############Parametros Epidemicos#####################"
   write(*,*) "Probabilidade de sitio estar infectado ", p0, "."
   write(*,*) ""
   write(*,*) "Taxa de recuperacao ", mu, "."	
   write(*,*) ""
   write(*,*) "Taxa de enfraquecimento imunologico ", alp, "."
   write(*,*) ""
   write(*,*) "Taxa de infeccao inicial ", lamb0, "."
   write(*,*) ""
   write(*,*) "Resolucao da taxa de infeccao ", dlamb, "."
   write(*,*) ""
   write(*,*) "Taxa de infeccao final ", lambdaf, "." 
   write(*,*) ""
   write(*,*) "Numero de pontos de taxa infeccao final ", nlamb, "."
   write(*,*) ""
   write(*,*) "#########################################################"
!#######################################################################
   call calcula_P_grau(rede)
!#######################################################################   
   !====================================================================
   if( nargus == 10)then
      nomeArquivo = trim(adjustl(local))//'lbd_vs_Y4_'//trim(adjustl(teoriaCM))//'_lamb_index_0'
   elseif( nargus == 9)then
      nomeArquivo = trim(adjustl(local))//'lbd_vs_Y4_'//trim(adjustl(teoriaCM))
   endif
   !====================================================================
   ! Testo se os arquivos acima existem.
   ! Logo abaixo, faco leituras, caso eles existam.
   ! Em caso negativo, eu chamo as rotinas apropriadas.
   !====================================================================
   inquire( file = trim(adjustl(nomeArquivo))//'.dat', exist = existe10)
   !====================================================================
   if( .not. existe10 )then
      write(*,*) "Arquivo nao existe"
      call aloca_listas_e_matrizes(rede)
      call metodo_bissecao(rede, alp, mu, lamb0, tole)
   else
      open(unit = 71, file = trim(adjustl(nomeArquivo))//'.dat', status = 'old') 
      read(71, *, iostat = st10) lamb0, Y4
      if( st10 /= 0)then
         write(*,*) "Nao houve leitura"
         call aloca_listas_e_matrizes(rede)
         call metodo_bissecao(rede, alp, mu, lamb0, tole)
      endif
      close(71)
   endif
   !====================================================================
   if( ( .not. existe10 ) .or. ( st10 /= 0))then
      !=================================================================   
      if( nargus == 10)then
         !==============================================================
         nomeArquivo = trim(adjustl(local))//'lbd_vs_Y4_'//trim(adjustl(teoriaCM))//'_lamb_index_0'
         open(333, file=trim(adjustl(nomeArquivo))//'.dat', status = 'unknown')
         !==============================================================
         nomeArquivo = trim(adjustl(local))//'lbd_vs_rho_'//trim(adjustl(teoriaCM))//'_lamb_index_0'
         open(334, file=trim(adjustl(nomeArquivo))//'.dat', status = 'unknown') 
         !==============================================================
      elseif( nargus == 9)then
         !==============================================================
         nomeArquivo = trim(adjustl(local))//'lbd_vs_Y4_'//trim(adjustl(teoriaCM))
         open(333, file=trim(adjustl(nomeArquivo))//'.dat', access = 'append', status = 'unknown')
         !==============================================================
         nomeArquivo = trim(adjustl(local))//'lbd_vs_rho_'//trim(adjustl(teoriaCM))
         open(334, file=trim(adjustl(nomeArquivo))//'.dat', access = 'append', status = 'unknown')
         !==============================================================
      endif
      !=================================================================
      write(*,*) "Ou nao existe ou leu ruim"
      !=================================================================
      write(333,*) lamb0, Y4
      close(333)
      !=================================================================
      write(334,*) lamb0, 0.0d0
      close(334)
      !=================================================================
      if( (nargus == 10) .and. (ind_lamb == 0) ) stop
      !=================================================================
   endif
   !====================================================================
   ! A principio
   !====================================================================
   if( nargus == 9 )then
      !=================================================================
      inquire(unit = 333, opened = taAberta)
      if(taAberta) close(333)   
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//trim(adjustl(teoriaCM))         
      !=================================================================
      open(333, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')   
      !=================================================================
      inquire(unit = 334, opened = taAberta)
      if(taAberta) close(334)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))
      !=================================================================
      inquire( file = trim(adjustl(nomeArquivo))//'.dat', exist = existe10)
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !=================================================================
      lambdaf = 0.0d0
      !=================================================================
      lerho:do
         read(334, *, iostat = st) lamb, rho
         if(st /= 0) exit lerho
         if( lamb > lambdaf) lambdaf = lamb
      enddo lerho
      if(lambdaf == 0.0d0)then
         lambdaf = lamb0 + dlamb * 500.0d0
         lamb = lamb0
      else
         lamb = lambdaf + dlamb
         lambdaf = lamb + dlamb * 500.0d0
      endif
      !=================================================================
      close(334)
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')
      !=================================================================
      inquire(unit = 335, opened = taAberta)
      if(taAberta) close(335)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//trim(adjustl(teoriaCM))         
      !=================================================================
      open(335, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')
      !=================================================================
      !-----------------------------------------------------------------
   elseif( nargus == 10 )then
      !-----------------------------------------------------------------
      inquire(unit = 333, opened = taAberta)
      if(taAberta) close(333)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
      !=================================================================
      open(333, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
      !=================================================================
      inquire(unit = 334, opened = taAberta)
      if(taAberta) close(334)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))       
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !=================================================================
      inquire(unit = 335, opened = taAberta)
      if(taAberta) close(335)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb)) 
      !=================================================================
      open(335, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !================================================================= 
      lambdaf = lamb0 + 1.0d0 * ind_lamb * dlamb
      lamb = lambdaf
   endif
   !====================================================================
   ! Eu testo se os arquivos supracitados existem ou contem
   ! dados apropriados. Em caso negativo, chamo os modulos
   ! apropriados e agora eu os escrevo.
   !====================================================================
   if( rede%nodes == 3*10**4 )then
      
      if(allocated(fi))then
         call calcula_distancias_ao_hub(rede)
         write(*,*) "Chamou calcula_distancias_ao_hub"      
         nomeArquivo = 'fi_pQMF_gm_'//trim(adjustl(gama_char))//'_tam_'//trim(adjustl(tam_char))//'_alp_'//trim(adjustl(alp_char2))//'_ams_'//trim(adjustl(indice))
         open(999, file = trim(adjustl(nomeArquivo))//'.dat', status = 'unknown')
         do i1 = 1, size(fi)
            write(999,*) rede%deg(i1), lista_distancias(i1), (fi(i1))**2.0d0
         enddo
         close(999)
         call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
         call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
      
         deallocate(lista_distancias)
      endif
   endif
   !====================================================================
   call aloca_listas(rede)
   !====================================================================
   teoriaCM = trim(adjustl(teoriaCM))
   usouCopia = .False.

   !stop "Calculou o limiar. Nao vai rodar dinamica"

   if(allocated(Ii_tild)) deallocate(Ii_tild)
   allocate(Ii_tild(rede%nodes))

   if(allocated(Ri_tild)) deallocate(Ri_tild)
   allocate(Ri_tild(rede%nodes))

   if(allocated(RIij_tild)) deallocate(RIij_tild)
   allocate(RIij_tild(rede%sumDeg))

   if(allocated(RRij_tild)) deallocate(RRij_tild)
   allocate(RRij_tild(rede%sumDeg))

   if(allocated(SIij_tild)) deallocate(SIij_tild)
   allocate(SIij_tild(rede%sumDeg))

llbd: do while( lamb <= lambdaf )
      !=================================================================
      ! tempo
      !=================================================================
      write(lamb_char, '(f10.6)') lamb     
      !=================================================================
      n_it_rho = 0
      n_it_Y4 = 0
      n_sub = 0
      !=================================================================
      if( nargus == 10)then
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
         !=================================================================
         open(336, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
         !=================================================================
         open(337, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rec_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
         !=================================================================
         open(339, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
      elseif( nargus == 9)then
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_'))//trim(adjustl(teoriaCM))
         !=================================================================
         open(336, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_'))//trim(adjustl(teoriaCM))
         !=================================================================
         open(337, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================         
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rec_'))//trim(adjustl(teoriaCM))
         !=================================================================
         open(339, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
         !=================================================================
      endif
      !=================================================================
      call checaBackUp(rede, trim(adjustl(local)))
      !=================================================================
      sumIi2 = sum(Ii)
      if( sumIi2 == 0.0d0)then
         call condicao_inicial(rede, p0, t0)
         rho0 = Ii(1)
         rho = rho0
         reco = Ri(1)
         t = 0.0d0
         Y4_old = 1.0_dp/(1.0_dp * comp_gigante)
         Y4 = Y4_old
         tempoEscrita = tempoEscrita0
         write(*,*) "Iniciou condicoes iniciais, nao backup"
      else
         reco = sum(Ri)/comp_gigante
      endif
      !=================================================================
lt:   do while(t <= tf)
         !==============================================================
         call tempoDeEscrever()
         !==============================================================
         write(336, *) t, rho
         write(337, *) t, Y4
         write(339, *) t, reco
         !==============================================================
         ! Aqui decidimos se o dt vai ou nao estragar tudo.
         ! Iniciamos Ii_tilda com um valor > 1.0d0
         !==============================================================
         choose_dt:do i1 = 1, size(vec_dt)
            dt = vec_dt(i1)
            call rk4_teste(dt, t, rede, alp, lamb, mu)
            !===========================================================
            rho = 0.0d0; reco = 0.0d0
            Ii_neg = .False.
            Ii_gt_1 = .False.
            !===========================================================
            conta_Ii:do i2 = 1, size(Ii_tild)   ! Calcula rho
               if(lista_de_clusters(i2) /= i_comp_gigante) cycle
               if( (Ii_tild(i2) < 0.0d0) .or. (Ri_tild(i2) < 0.0d0) .or. ( (1.0d0 - Ii_tild(i2) - Ri_tild(i2)) <= 0.0d0 ) )then
                  Ii_neg = .True.
                  exit conta_Ii
               elseif( (Ii_tild(i2) > 1.0d0) .or. (Ri_tild(i2) > 1.0d0) .or. ( ( 1.0d0 - Ii_tild(i2) - Ri_tild(i2)) > 1.0d0 ))then
                  Ii_gt_1 = .True.
                  exit conta_Ii
               endif
               rho = rho + Ii_tild(i2)
               reco = reco + Ri_tild(i2)
            enddo conta_Ii
            !===========================================================
            if(.not. Ii_neg  )then
               if( .not. Ii_gt_1 )then
                  Ii = Ii_tild
                  Ri = Ri_tild
                  RIij = RIij_tild
                  RRij = RRij_tild
                  SIij = SIij_tild
                  exit choose_dt
               endif
            endif
         enddo choose_dt
         !if(i1 > 1) write(*,*) "dt atualizado para = ", dt, "de indice = ", i1
         !==============================================================
         !call rk4(dt, t, rede, rede%nodes, rede%sumDeg, alp, lamb, mu)
         !==============================================================
  
         if( ( Ii_neg == .True. ) .or. ( Ii_gt_1 == .True. ) ) stop "Esgotou a precisao da Maquina, sem obter resultados fisicos"
         !==============================================================
         rho  = rho/( 1.0_dp * comp_gigante )
         reco  = reco/( 1.0_dp * comp_gigante )         
         
         sumIi2 = (sum(Ii**2.0_dp))**0.5_dp
         
		 Y4 = sum( (Ii/sumIi2)**4.0_dp )
         !##############################################################
         t = t + dt
         !##############################################################
         
         if( abs(rho - rho0)/rho0 < tole )then
            n_it_rho = n_it_rho + 1
            if( n_it_rho >= int(5.0d0/dt0) )then
               write(*,*) "Conv. lbd = ", lamb, " e rho = ", rho, "."
               exit lt
            endif
         elseif( abs(Y4 - Y4_old)/Y4_old < tole )then
            n_it_Y4 = n_it_Y4 + 1
            if( n_it_Y4 >= int(5.0d0/dt0) )then
               write(*,*) "Conv. lbd = ", lamb, " e Y4 = ", Y4, "."
               exit lt
            endif
         else
            n_it_rho = 0; n_it_Y4 = 0
         endif
         !==============================================================
         !if( t > 5.0d0)then
            !do i1 = 1, rede%nodes
               !do i12 = rede%aux(i1), rede%aux(i1) + rede%deg(i1) - 1
                  !i2 = rede%listAdj(i12)
                  !i21 = spec(i12)
                  !write(*,*) "Ii - SIij - RIij - IIij = ", abs(Ii(i1) - Ii(i2) - SIij(i21) - RIij(i21)  + SIij(i12) + RIij(i12))
               !enddo            
            !enddo
            !stop
         !endif
         !==============================================================
         !   Se deu ruim por causa de dt grande, a gente pega uma
         !   configuracao passada e abaixa o dt
         !==============================================================
         !   Antes faziamos assim :/ ...
         !   Isso devia ir pra um museu!
         !==============================================================
         !stop "Ajuste dt para que rho deixe de ser negativo."
         !if( rho > 1.0_dp) stop "Ajuste dt para que rho deixe de ser > 1.0"
         !##############################################################
         rho0 = rho 
         Y4_old = Y4        
      enddo lt
      !=================================================================
      close(336)
      close(337)
      !=================================================================
      write(*,*) "lambda e rho = ", lamb, rho      
      !=================================================================
      if( (n_it_Y4 >= int(5.0d0/dt) ) .or. (n_it_rho >= int(5.0d0/dt) ) )then
         write(333, *) lamb, Y4
         write(334,*) lamb, rho   
         write(335,*) lamb, t                     
      endif
      !=================================================================
      lamb = lamb + dlamb
      !=================================================================
      if(nargus == 9)then
         call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_global.dat')
         call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_global.tar.gz')
      elseif(nargus == 10)then
         call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
         call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.tar.gz')
      endif
      !=================================================================
   enddo llbd
   close(333)
   close(334)
   close(335)
   
!#######################################################################   

      if(allocated(k1_Ii))then
         deallocate(k1_Ii) 
         deallocate(k2_Ii)
         deallocate(k3_Ii)
         deallocate(k4_Ii)

         deallocate(k1_Ri) 
         deallocate(k2_Ri)
         deallocate(k3_Ri)
         deallocate(k4_Ri)

         deallocate(k1_SIij) 
         deallocate(k2_SIij)
         deallocate(k3_SIij)
         deallocate(k4_SIij)

         deallocate(k1_RIij) 
         deallocate(k2_RIij)
         deallocate(k3_RIij)
         deallocate(k4_RIij)

         deallocate(k1_RRij) 
         deallocate(k2_RRij)
         deallocate(k3_RRij)
         deallocate(k4_RRij)
      endif            

  contains

   function Delta_Ii_Euler_hub(this, ihub, n_sitios, n_arestas, dt, alf, lbd, mu, Ii, SIij)
      class(grafo), intent(in) :: this
      integer, intent(in) :: ihub
      integer, intent(in) :: n_sitios
      integer (kind=8), intent(in) :: n_arestas
      real(dp) :: Delta_Ii_Euler_hub
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: alf
      real(dp), intent(in) :: lbd
      real(dp), intent(in) :: mu
      real(dp) :: Ii(n_sitios)
      real(dp) :: SIij(n_arestas)
      integer :: k1
      integer(kind = 8) :: k12
      real(dp) :: SIij_aux

      !#################################################################
      k1 = ihub
      !#################################################################
         SIij_aux = 0.0_dp
         !##############################################################
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            SIij_aux = SIij_aux + SIij(k12)
         enddo
         !##############################################################
         Delta_Ii_Euler_hub = dt * (lbd * SIij_aux - mu * Ii(ihub))
      !#################################################################
   end function   
   !====================================================================
      subroutine criaRedeEClassificaClusters(this, tam_rede1, grau_min1, grau_max1, gama_exp1, seme1)
        type(grafo_PL_UCM) :: this
        integer :: tam_rede1
        integer :: grau_min1
        real(dp) :: grau_max1
        real(dp) :: gama_exp1
        integer :: seme1
!#######################################################################   
        call this%iniciaGrafo(tam_rede1)
!#######################################################################
        call this%inicia(grau_min1, grau_max1, gama_exp1, seme1)
!#######################################################################
        call this%liga(seme1, .False.) 
!#######################################################################        
        call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat') 
!#######################################################################
      end subroutine
!#######################################################################

      subroutine entradaArgumentos()
         nargus = iargc()

         if(nargus == 9)then
            !#############################
            ! Amostra
            !#############################
            write(*,*) "Recebeu 9 arqumentos"
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_amostra

            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede

            !#############################
            ! Tamanho	da rede
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_min

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) gama_exp

            !#############################
            ! Lambda0
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor


            !#############################
            ! Lambdaf
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lambdaf

            !#############################
            ! Alfa
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
            
            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) int_soCalculaIPR
      
            if(int_soCalculaIPR == 0)then
               soCalculaIPR = .False.
            elseif(int_soCalculaIPR == 1)then
               soCalculaIPR = .True.
            else
               stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
            endif     
                        
         elseif(nargus == 10)then
            write(*,*) "Recebeu 9 arqumentos"
            !#############################
            ! Amostra
            !#############################
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_amostra

            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede

            !#############################
            ! Tamanho	da rede
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_min

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) gama_exp

            !#############################
            ! Lambda0
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor


            !#############################
            ! Lambdaf
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lambdaf

            !#############################
            ! Alfa
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
      
            !#############################
            ! Indice do ponto
            ! Quando fizermos
            ! paralelizacao burra.
            !#############################
            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_lamb   
            write(char_ind_lamb, '(I0)') ind_lamb
        
            char_ind_lamb = trim(adjustl(char_ind_lamb)) 

            call getarg(10, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) int_soCalculaIPR
      
            if(int_soCalculaIPR == 0)then
               soCalculaIPR = .False.
            elseif(int_soCalculaIPR == 1)then
               soCalculaIPR = .True.
            else
               stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
            endif
         else
            stop "Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
         endif

!#######################################################################
       	 if( trim(adjustl(resultados)) == 'Rst_'//trim(adjustl(teoriaCM))//'_Corte_Rigido' )then 
            call acha_cutoff_rigido(grau_min, gama_exp, tam_rede)
       	 elseif(trim(adjustl(resultados)) =='Rst_'//trim(adjustl(teoriaCM))//'_Corte_sqrtN' )then
            grau_max = (1.0_dp * tam_rede)**(0.5_dp)
         elseif(trim(adjustl(resultados)) == 'Rst_'//trim(adjustl(teoriaCM))//'_Corte_2sqrtN' )then
            grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
       	 endif
         dlamb = 0.0125_dp/(1.0_dp * divisor)
!#######################################################################
         p0 = 10.0_dp/(1.0_dp * tam_rede) !1.0d-2
!#######################################################################
!  Indice da semente
!#######################################################################
         write(gama_char,'(f5.2)') gama_exp
!#######################################################################
         write(indice,'(I0)') ind_amostra
!#######################################################################   
         if(tam_rede == 10)then
            tam_char = '10'
         elseif(tam_rede == 100)then
            tam_char = '100'
         elseif(tam_rede == 1000)then
            tam_char = '1k'
         elseif(tam_rede == 3000)then
            tam_char = '3k'
         elseif(tam_rede == 10000)then
            tam_char = '10k'
         elseif(tam_rede == 30000)then
            tam_char = '30k'
         elseif(tam_rede == 100000)then
            tam_char = '100k'
         elseif(tam_rede == 300000)then
            tam_char = '300k'
         elseif(tam_rede == 1000000)then
            tam_char = '1M'
         elseif(tam_rede == 3000000)then
            tam_char = '3M'
         elseif(tam_rede == 10000000)then
            tam_char = '10M'
         elseif(tam_rede == 30000000)then
            tam_char = '30M'
         elseif(tam_rede == 100000000)then
            tam_char = '100M'
         else
            stop 'Escolha um tamanho de rede dentro do catalogo'
         endif         
      end subroutine

!=======================================================================
   subroutine calcula_P_grau(this)      
      class(grafo) :: this
      integer :: i1
!#######################################################################
      if(allocated(P_grau)) deallocate(P_grau)
      allocate(P_grau(this%degMin:this%degMax))
      P_grau = 0_dp
      do i1 = 1, this%nodes
         P_grau(this%deg(i1)) = P_grau(this%deg(i1)) + 1.0_dp
      enddo
      P_grau = P_grau/(1.0_dp * this%nodes)
      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_P_grau_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      open(800, file=trim(adjustl(arq_1)), status='unknown')
      do i1 = this%degMin, this%degMax
         if(P_grau(i1) == 0.0_dp) cycle
         write(800,*) i1, P_grau(i1)
      enddo
      close(800)
      deallocate(P_grau)
!#######################################################################      
   end subroutine
!=======================================================================
   subroutine kNN_e_clustering(this)
      class(grafo) :: this

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call calcula_k_nn(this, .True., 800, trim(adjustl(arq_1)))
      close(800)

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call clustering(this,.True., 800, trim(adjustl(arq_1)))
      close(800)
   end subroutine
!=======================================================================   
	subroutine checaBackUp(this, pasta_fonte)
		 !--------------------------------------------------------------
		 class(grafo) :: this
		 character(len = *) :: pasta_fonte
		 character(len = 1000) :: nomeArquivo
		 integer :: i1, i2, i3
		 integer(kind = 8) :: i12, i21, i23, i13, i31, i32
		 logical :: existe
		 integer :: st
		 real(dp), allocatable :: rho_i1(:)
		 real(dp) :: norm_rhoi1
		 character(len = 100) :: C_lamb_f
		 character(len = 100) :: C_tempo_f, C_tempoEscrita_f
		 character(len = 100) :: C_Ii_f, C_Ri_f, C_SIij_f, C_RIij_f, C_RRij_f
		 real(dp) :: t_lido, tempoEscrita_lido
		 !--------------------------------------------------------------
		 pasta_fonte = trim(adjustl(pasta_fonte))
		 !--------------------------------------------------------------
		 if(nargus == 10 )then
			C_lamb_f = trim(adjustl(C_lamb))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_tempo_f = trim(adjustl(C_tempo))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_tempoEscrita_f = trim(adjustl(C_tempoEscrita))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_Ii_f = trim(adjustl(C_I_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_Ri_f = trim(adjustl(C_R_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_SIij_f = trim(adjustl(C_SIij))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_RIij_f = trim(adjustl(C_RIij))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_RRij_f = trim(adjustl(C_RRij))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
		 elseif( nargus == 9)then
			C_lamb_f = trim(adjustl(C_lamb))//'_global'
			C_tempo_f = trim(adjustl(C_tempo))//'_global'
			C_tempoEscrita_f = trim(adjustl(C_tempoEscrita))//'_global'
			C_Ii_f = trim(adjustl(C_I_i))//'_global'
			C_Ri_f = trim(adjustl(C_R_i))//'_global'
			C_SIij_f = trim(adjustl(C_SIij))//'_global'
			C_RIij_f = trim(adjustl(C_RIij))//'_global'
			C_RRij_f = trim(adjustl(C_RRij))//'_global'
		 endif
		 Ii = 0.0d0
		!C_lamb_f, C_tempo_f, C_Ii_f, C_Ri_f, C_SIij_f, C_RIij_f, C_RRij_f
		!===========================================================
		!-------------------Leitura de LAMBDA-----------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_lamb_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe ) 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

		open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
		read(776,*, iostat=st) lbd_lido
		close(776)
		if( st /= 0 )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		!===========================================================
		!-------------------Leitura do TEMPO------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_tempo_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe )
        !===============================================================
        ! C_tempo_f nao existe
        !===============================================================

		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		
		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

		open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
		read(776,*, iostat=st) t_lido
		close(776)
		if( st /= 0 )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		!===============================================================
		!-------------------Leitura do TEMPO DE ESCRITA-----------------
		!===============================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_tempoEscrita_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe )
        !===============================================================
        ! C_tempo_f nao existe
        !===============================================================

		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		
		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

		open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
		read(776,*, iostat=st) tempoEscrita_lido
		close(776)
		if( st /= 0 )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		!===========================================================
		!-------------------Leitura de I_i--------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_Ii_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe ) 

		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)
        
        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   !============================================================
		   ! Eu zeraria Ii aqui, mas ele ainda nao foi lido.
		   ! Portanto, ele jah eh zero.
		   !============================================================
		   return
        endif
        !===============================================================
		do i1 = 1, size(Ii)
		   read(776,*, iostat=st) Ii(i1)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !============================================================
		      ! Eh preciso zerar Ii, afim de que meu programa entenda
		      ! que nao foi feita uma leitura apropriada do backup
		      ! e ele chame as condicoes iniciais apropriadas
		      !============================================================
			  Ii = 0.0d0
			  return
		   endif
		enddo
		close(776)
		!===========================================================
		!-------------------Leitura de R_i--------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_Ri_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe ) 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')
		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)

        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
        endif
        !===============================================================
		do i1 = 1, size(Ri)
		   read(776,*, iostat=st) Ri(i1)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !=========================================================
		      ! Preciso zerar ao menos Ii, afim de que meu programa
		      ! entenda que checaBackup nao achou backup
		      ! e chame a rotina condicao_inicial.
		      !=========================================================
		      Ii = 0.0d0
		      return
		   endif
		enddo
		close(776)
		!===============================================================
		!-------------------Leitura de SI_ij-------------------------
		!===============================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_SIij_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
		 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   !============================================================
		   ! Preciso zerar pelo menos o Ii, 
		   ! senao, o programa vai entender que houve leitura correta.
		   !============================================================
		   Ii = 0.0d0
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)
        
        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
        endif
        !===============================================================
		do i12 = 1, size(SIij)
		   read(776,*, iostat=st) SIij(i12)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !=========================================================
		      ! Preciso zerar ao menos Ii, pois foi a condicao
		      ! que eu escolhi para declarar que checaBackup
		      ! nao teve sucesso.
		      !=========================================================
			  Ii = 0.0d0
			  return
		   endif
		enddo
		close(776)
		!===============================================================
		!-------------------Leitura de RI_ij-------------------------
		!===============================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_RIij_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
		 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   !============================================================
		   ! Preciso zerar pelo menos o Ii, 
		   ! senao, o programa vai entender que houve leitura correta.
		   !============================================================
		   Ii = 0.0d0
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)
        
        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
        endif
        !===============================================================
		do i12 = 1, size(RIij)
		   read(776,*, iostat=st) RIij(i12)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !=========================================================
		      ! Preciso zerar ao menos Ii, pois foi a condicao
		      ! que eu escolhi para declarar que checaBackup
		      ! nao teve sucesso.
		      !=========================================================
			  Ii = 0.0d0
			  return
		   endif
		enddo
		close(776)
        !===============================================================
        ! !C_lamb_f, C_tempo_f, C_Ii_f, C_Ri_f, C_SIij_f, C_RIij_f, C_RRij_f
        !===============================================================
        
		!===============================================================
		!-------------------Leitura de RR_ij-------------------------
		!===============================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_RRij_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
		 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   !============================================================
		   ! Preciso zerar pelo menos o Ii, 
		   ! senao, o programa vai entender que houve leitura correta.
		   !============================================================
		   Ii = 0.0d0
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)
        
        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
        endif
        !===============================================================
		do i12 = 1, size(RRij)
		   read(776,*, iostat=st) RRij(i12)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !=========================================================
		      ! Preciso zerar ao menos Ii, pois foi a condicao
		      ! que eu escolhi para declarar que checaBackup
		      ! nao teve sucesso.
		      !=========================================================
			  Ii = 0.0d0
			  return
		   endif
		enddo
		close(776)
		!===============================================================
        ! Agora que leu tudo, sera feita a estatistica!
		!===============================================================
		rho = 0.0d0
		do i1 = 1, size(Ii)
		   if( lista_de_clusters(i1) /= i_comp_gigante) cycle
		   rho = rho + Ii(i1)
		enddo
		rho = rho/(1.0d0 * comp_gigante)
		!---------------------------------------------------------------    
		!if(allocated(rho_i1)) deallocate(rho_i1)
		!allocate(rho_i1(size(Ii)))
		!---------------------------------------------------------------
		! Calcular o NAV via mensagens que chegam dava um resultado
		! ligeiramente diferente para o Y4. Levemente maior.
		! Isso acontecia porque eu nao estava levando em conta a
		! equacao
		! \rho_i = (lambda/mu) * s_i * \sum_j A_{ij} * I_{ji}.
		! Apos fazer isso, a concordancia ficou muito melhor.
		!---------------------------------------------------------------
		!rho_i1 = 0.0d0
		!---------------------------------------------------------------
		!do i1 = 1, size(Ii)
		!   if( lista_de_clusters(i1) /= i_comp_gigante) cycle
		!   do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
		!	  i21 = link_back(i12)
		!	  rho_i1(i1) = rho_i1(i1) + Iji(i21)
		!   enddo
		!   rho_i1(i1) = lbd_lido * (1.0d0 - Ii(i1) - Ri(i1) )* rho_i1(i1)
		!enddo
		!norm_rhoi1 = 0.0d0

		!norm_rhoi1 = (sum(rho_i1**2.d0))**0.5d0
		!rho_i1 = rho_i1/norm_rhoi1
		!Y4 = sum(rho_i1**4.0d0)
		
		!write(*,*) t_lido, Y4

		!===============================================================
		!norm_rhoi1 = (sum(Ii**2.0d0))**0.5d0
		
		!Y4 = sum((Ii/norm_rhoi1)**4.0d0)
		!write(*,*) t_lido, Y4
		!===============================================================
		!deallocate(rho_i1)
		!stop
		!===============================================================
		norm_rhoi1 = (sum(Ii**2.0d0))**0.5d0
		
		Y4 = sum((Ii/norm_rhoi1)**4.0d0)
		!===============================================================
		Y4_old = Y4
		rho0 = rho
		!===============================================================
		t = t_lido
		lamb = lbd_lido
	end subroutine
!=======================================================================
   subroutine TempoDeEscrever()
                         
         if( t >= tempoEscrita)then
            !===========================================================       
            do i1 = 1, size(Ii)
               if( (Ii(i1) > 1.0d0) .or. (Ii(i1) < 0.0d0) ) return
            enddo
            !===========================================================
            do i1 = 1, size(Ri)
               if( ( Ri(i1) > 1.0d0) .or. (Ri(i1) < 0.0d0) )return
            enddo            
            !===========================================================
            do i12 = 1, size(SIij)
               if( ( SIij(i12) > 1.0d0) .or. (SIij(i12) < 0.0d0) )return
            enddo
            !===========================================================
            do i12 = 1, size(RIij)
               if( ( RIij(i12) > 1.0d0) .or. (RIij(i12) < 0.0d0) )return
            enddo
            !===========================================================
            do i12 = 1, size(RRij)
               if( ( RRij(i12) > 1.0d0) .or. (RRij(i12) < 0.0d0) )return
            enddo                              
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_i))//'_global'
            endif
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')            
            !===========================================================            
            do i1 = 1, size(Ii)
               write(777,*) Ii(i1)
            enddo
            !===========================================================
            close(777)                 
            !===========================================================
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================

            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_R_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_R_i))//'_global'
            endif            
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================            
            do i1 = 1, size(Ri)
               write(777,*) Ri(i1)
            enddo
            !===========================================================
            close(777)
            !===========================================================
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            ! SI_ij
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_SIij))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_SIij))//'_global'
            endif
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================               
            do i12 = 1, size(SIij)
               write(777,*) SIij(i12)
            enddo
            !===========================================================
            close(777)
            !===========================================================            
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            ! RI_ij
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_RIij))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_RIij))//'_global'
            endif
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================               
            do i12 = 1, size(RIij)
               write(777,*) RIij(i12)
            enddo
            !===========================================================
            close(777)
            !===========================================================            
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            ! RR_ij
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_RRij))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_RRij))//'_global'
            endif
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================               
            do i12 = 1, size(RRij)
               write(777,*) RRij(i12)
            enddo
            !===========================================================
            close(777)
            !===========================================================            
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            ! Tempo
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempo))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempo))//'_global'
            endif           
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) t
            close(777)
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_lamb))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_lamb))//'_global'
            endif  
            !===========================================================       
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) lamb
            close(777)
            !===========================================================
            tempoEscrita = tempoEscrita + tempoEscrita0
            !===========================================================
            ! Tempo de Escrita
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempoEscrita))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempoEscrita))//'_global'
            endif           
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) tempoEscrita
            close(777)
            !===========================================================
         endif
   end subroutine
!=======================================================================
   function g_func(k_s, gam1, Nstr)
      real(dp) :: k_s
      real(dp) :: gam1
      integer :: Nstr
      real(dp) :: g_func
      real(dp) :: Ap
      !#################################################################
      !   Quando g_func for usada pela primeira vez,
      !   gbuffer = 0.0d0 e kmin = this%degMin,
      !   k_s recebe o valor que se quer aproximar de kc
      !   
      !#################################################################      
      !g_func = gbuffer
                  
      Ap = Ap_list(int(k_s))
      
      Ap = 1.0_dp/Ap
      
      g_func = k_s - Ap**(1.0_dp/gam1) * (1.0_dp * Nstr)**(1.0_dp/gam1)
      
   end function                          

   subroutine acha_cutoff_rigido(kming, gam_p, N_p)
      integer :: kming
      real(dp) :: gam_p
      integer :: N_p
      real(dp), parameter :: tol = 5d-5
      real(dp) :: gminus, gmais
      real(dp) :: kminus, kmais
      real(dp) :: k_p
      real(dp) :: gp
      real(dp) :: Ap1
      integer :: kl1, kl2
      integer, parameter :: N_it = 10**3
 
      !#################################################################
      !   Inicio
      !#################################################################           
      kminus = 1.0_dp * kming
      kmais = 1.5_dp * kming * (1.0_dp * N_p)**(1.0_dp/gam_p)
      
      if(allocated(Ap_list)) deallocate(Ap_list)
      allocate(Ap_list(int(kming):(int(kmais))))
      
      gminus = g_func(kminus, gam_p, N_p)
      
      if(gminus >= 0.0_dp) stop "Precisamos de um valor de kminus para que gminus < 0"
       
      Ap1 = 0.0_dp
      do kl1 = kming, int(kmais)
         Ap1 = Ap1 + (1.0_dp * kl1)**(-gam_p)
         Ap_list(kl1) = Ap1
      enddo

      gmais = g_func(kmais, gam_p, N_p)
   
      if(gmais <= 0.0_dp) stop "Precisamos de um valor para kmais, tal que gmais > 0"
      !#################################################################
      !   Execucao
      !#################################################################
      kl1 = 1
      do while(kl1 <= N_it)
         k_p = kminus + (kmais - kminus)/2.0_dp
         gp = g_func(k_p, gam_p, N_p)
         if((gp == 0.0_dp).or.((kmais-kminus)/2.0_dp <= tol))then
            grau_max = k_p
            write(*,*) "Achou o grau max apos N = ", kl1, " iteracoes."
            exit
         endif
         kl1 = kl1 + 1
         if(gminus * gp > 0.0_dp)then
            kminus = k_p
            gminus = gp
         else
            kmais = k_p
         endif
      enddo      
      !#################################################################
      !   Final
      !#################################################################
      deallocate(Ap_list)                   
   end subroutine
   
  
  
end program
