module limiaresCampoMedio
   
   use geraRede
   use mod_rndgen
   
   implicit none
   
   integer, private, parameter :: dp = kind(0.0d0)
   
   real(dp), allocatable :: x(:), y(:)
   real(dp) :: Y4
   real(dp), allocatable :: f1(:)   
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

   lambda_mean = dlambda0
   
   lambda_minus = 0.0d0
   
   lambda_plus = 0.0d0
   !====================================================================
   do while(.True.)
   
      i1 = i1 + 1
      write(*,*) "Entrou no metodo da potencia"
!=======================================================================
     !call metodo_potencia(this, alp1, mu1, lambda1, autovalor1, tole)
      call metodo_potencia(this, alp, mu, lambda_mean, autovalor, tole)
!=======================================================================
      !write(*,*) "Saiu do metodo da potencia com autovalor = ", autovalor
      !write(789, *) i1, lambda_mean

      if( ((lambda_minus * lambda_plus) == 0.0d0) )then
         !--------------------------------------------------------------
         !write(*,*) "Lambda minus * lambda plus == 0"
         !write(*,*) "Lambda minus = ", lambda_minus
         !write(*,*) "lambda_plus = ", lambda_plus
         !--------------------------------------------------------------         
         if( autovalor < 0.0d0)then
            lambda_minus = lambda_mean
            autovalor_minus = autovalor
            lambda_mean = lambda_mean + dlambda
         else
            !write(*,*) "Entrou no Else"            
            if( autovalor < tole )then
               write(*,*) 'Lambda_mean convergiu para = ', lambda_mean
               close(789)
               exit
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
            exit
         else
            lambda_minus = lambda_mean
            autovalor_minus = autovalor
         endif
         lambda_mean = (lambda_minus + lambda_plus)/2.0d0
!=======================================================================
         if( abs( lambda_plus - lambda_minus ) < tole)then
            write(*,*) 'Lambda_mean convergiu para = ', lambda_mean
            close(789)
            exit
         endif
!=======================================================================      
      endif
!=======================================================================      
   enddo
!=======================================================================

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
         class is (grafo)
            write(*,*) "Deu classe grafo"            
            qm = 1.0d0 * this%degMean               
            limiarQMF = 1.0d0/(1.0d0 * qm) 
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

      if(allocated(f1)) deallocate(f1)
      
      allocate(f1(this%nodes))            
 
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
            f1 = (y/y_norm + x/x_norm)/2.0d0
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
   real(dp), allocatable :: SIij(:), RRij(:), RIij(:)
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
      integer :: k1, k2, k12, k3, k13, k23
      
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
      
      if(allocated(v1)) deallocate(v1)
      allocate(v1(this%nodes))
      
   end subroutine
   
   subroutine condicao_inicial(this, I_0, t1)     

      class(grafo), intent(in) :: this
      real(dp), intent(in) :: I_0
      integer:: k1, k2, k12, k21, k23
      logical :: isI_i
      integer :: iost
      character(len=30) :: Nichar
      character(len=10) :: lambdachar
      character(len=100) :: arquivo_char2
      integer :: sum_deg
      real(dp) :: t1    
      t = t1
                     
      do k1 = 1, this%nodes
         if(lista_de_clusters(k1) /= i_comp_gigante)then
            Ii(k1) = 0.0_dp
            Ri(k1) = 0.0_dp
            cycle
         endif
         Ii(k1) = I_0
         v1(k1) = Ii(k1)
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
						
!#######################################################################
        !---------------------------------------------------------------
        Ii =   Ii   + (1.0_dp/6.0_dp) * (k1_Ii   + 2.0_dp * k2_Ii   + 2.0_dp * k3_Ii   + k4_Ii)
        !---------------------------------------------------------------
        v1 = Ii
        !---------------------------------------------------------------
        Ri =   Ri   + (1.0_dp/6.0_dp) * (k1_Ri   + 2.0_dp * k2_Ri   + 2.0_dp * k3_Ri   + k4_Ri)
        !---------------------------------------------------------------
        RIij = RIij + (1.0_dp/6.0_dp) * (k1_RIij + 2.0_dp * k2_RIij + 2.0_dp * k3_RIij + k4_RIij)
        !---------------------------------------------------------------        
        RRij = RRij + (1.0_dp/6.0_dp) * (k1_RRij + 2.0_dp * k2_RRij + 2.0_dp * k3_RRij + k4_RRij)
        !---------------------------------------------------------------        
        SIij = SIij + (1.0_dp/6.0_dp) * (k1_SIij + 2.0_dp * k2_SIij + 2.0_dp * k3_SIij + k4_SIij) 
        !---------------------------------------------------------------        
!#########################SUBROTINA##############################################################      
   end subroutine
   !####################################################################
   !   Funcoes
   !####################################################################

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
      integer :: k1, k12
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
      integer :: k1, k2, k3, k12, k21, k13, k31, k23, k32
      
      real(dp) :: SIij_aux
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
            SIij_aux = 0.0_dp
            !###########################################################
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               if( k23 == k21 ) cycle
               SIij_aux = SIij_aux + SIij(k23)
            enddo
            !###########################################################
            f_SIij(k12) = f_SIij(k12) + lbd * ( SSij_aux/ ( 1.0_dp - Ii(k2) - Ri(k2) ) ) * SIij_aux
            !###########################################################
            SIij_aux = 0.0_dp
            do k13 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
               if( k13 == k12 ) cycle
               SIij_aux = SIij_aux + SIij(k13)
            enddo
            !###########################################################
            ! Antes tinha um sinal de -. Optei por mudar o sinal de Si.
            f_SIij(k12) = f_SIij(k12) - lbd * ( SIij(k12) /( 1.0_dp - Ri(k1) - Ii(k1) ) ) * SIij_aux
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
      integer :: k1, k2, k12, k21
      
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
      integer :: k1, k2, k3, k12, k21, k13, k31, k23, k32
      
      real(dp) :: SIij_aux
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
            !write(*,*) abs( Ii(k2) - SIij(k12) - RIij(k12) - Ii(k1) + SIij(k21) + RIij(k21) )
            !###########################################################
            SIij_aux = 0.0_dp
            !###########################################################
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               if( k23 == k21 ) cycle               
               SIij_aux = SIij_aux + SIij(k23)               
            enddo
            !###########################################################
            f_RIij(k12) = f_RIij(k12) + lbd * ( ( Ri(k1) -RIij(k12) - RRij(k12) )/( 1.0_dp - Ii(k2) - Ri(k2) ) ) * SIij_aux
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
   type(grafoRRN_Plus_Star) :: rede
!#######################################################################   
   real(dp) :: dlamb
   real(dp) :: lamb0, lambdaf
   real(dp) :: lamb
   integer :: nlamb = 1000
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
   real(dp) :: rho, rho0
   real(dp) :: tol, tole
!#######################################################################
   !type(rndgen) :: gen1
   integer :: seed(10)
   integer :: seed1
   type(rndgen) :: ger_inic
!#######################################################################
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   integer :: grau_max      
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
   integer :: i1, i12, i2, i3, i4, i5, i6, i7, i8, i9 !, j1, j2, j3
   logical :: T_vs
!#######################################################################
   integer :: sumDeg2
   real(dp) :: t0, tf
   real(dp) :: dt_m
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   integer :: ntempo
   
   integer :: n_it, n_sub, per_conv
   
   real(dp) :: p0
   character(len=300) :: cwd, resultados, tipoCorte
   character(len=1000) :: local
   character(len=1000) :: nomeArquivo
   character(len=20) :: buffer
   !##############################################################
   
   real(dp) :: qm, q2m

   integer :: nargus, ind_lamb   
   character(len=3) :: char_ind_lamb
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
   integer :: grau_RRN
   integer :: j1
!#######################################################################   
   seed1=947361823
   salvouCentralidade = .False.
!#######################################################################   
 
   teoriaCM = trim(adjustl('pQMF'))
      
   resultados = 'Rst_'//trim(adjustl(teoriaCM))//'_RRN_Plus_Star'
    
   resultados = trim(adjustl(resultados))

   call system('mkdir -p '//resultados)
    
   local = trim(adjustl(resultados))//"/"

   call entradaArgumentos()

!#######################################################################
   !====================================================================
   !  Se der ruim no dt, ele eh dividido por dois.
   !====================================================================   
   dt = 1.0d-1 
   !====================================================================
   ! A principio, t0 = 0, mas, se houver um estado salvo, muda
   ! para o t salvo no arquivo.
   !====================================================================
   t0 = 0.0_dp
   tf = 10000.0_dp
   
   tempoEscrita0 = tf/100.0d0
   
   ntempo = int( (tf- t0)/dt )
   tole = 1d-7
   per_conv = int( 10.0_dp/dt )
!#######################################################################

!#######################################################################
   call ger_inic%init(seed1)
!   i2 = 1
!   do i1 = 1, 1000
!      if(mod(i1,100) > 0) cycle
!      seed(i2)  = ger_inic%int(100000000,999999999)
!      write(*,*) i1, seed(i2)
!      i2 = i2+1      
!   enddo

   i2 = 1
   do i1 = 1, 1000
      if(mod(i1,100) /= 50)then !Antes, eu dava um cycle. Agora eh assim      
         so_gasta = ger_inic%int(100000000,999999999)
      else
         seed(i2)  = ger_inic%int(100000000,999999999)
         write(*,*) i1, seed(i2)
         i2 = i2+1
      endif
   enddo

   local = trim(adjustl(trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
      
   call system('mkdir -p '//trim(adjustl(local)) )


!#######################################################################
!				Inicia grafo
!#######################################################################
 !  arquivo_rede_real='s00088.s838.net.edg'
 !  call rede%RedeReal(arquivo_rede_real, 111)
 !  close(111)
!#######################################################################
   write(*,*) "Grau RRN = ", grau_RRN
   write(*,*) "Grau_max=", grau_max
   
   call criaRedeEClassificaClusters(rede, tam_rede, grau_RRN, grau_max)

   !====================================================================
   write(*,*) "######################Dados da Rede######################"
   write(*,*) ""
   write(*,*) "Tamanho da rede ", rede%nodes, "."
   write(*,*) ""
   write(*,*) "Fracao correspondente aa componente gigante ", 100.0 * comp_gigante/rede%nodes,"%", "."
   write(*,*) ""
   write(*,*) "Grau minimo ", rede%degMin, ".", " Grau maximo ", rede%degMax, "."
   write(*,*) ""
!#######################################################################
   write(alp_char2, '(f9.3)') alp

   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//'/'))
 
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

!#######################################################################

!#######################################################################
   
   !====================================================================
   call aloca_listas(rede)
   !====================================================================
   usouCopia = .False.

   if( nargus == 9)then

      inquire(unit = 333, opened = taAberta)
      if(taAberta) close(333)   
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))
      !=================================================================
      !=================================================================
      open(333, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')   
      !=================================================================

      inquire(unit = 334, opened = taAberta)
      if(taAberta) close(334)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')
      !=================================================================  

      inquire(unit = 335, opened = taAberta)
      if(taAberta) close(335)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))         
      !=================================================================
      open(335, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')
      !=================================================================

   elseif( nargus == 10)then

      inquire(unit = 333, opened = taAberta)
      if(taAberta) close(333)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
      !=================================================================
      !=================================================================
      open(333, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
      !=================================================================
   
      inquire(unit = 334, opened = taAberta)
      if(taAberta) close(334)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !=================================================================

      inquire(unit = 335, opened = taAberta)
      if(taAberta) close(335)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))//'_lamb_index_'//trim(adjustl(char_ind_lamb)) 
      !=================================================================
      open(335, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !=================================================================             
   endif

   call aloca_listas_e_matrizes(rede)
   lamb0 = 1.0d0 * dlamb
   !====================================================================
   call metodo_bissecao(rede, alp, mu, lamb0, 1.0d-7)!TestaSeTemCopiaEstadoSistema(usouCopia)
   !====================================================================
   write(333,*) lamb0, Y4
   write(334,*) lamb0, 0.0d0
   
   if( rede%nodes == 10**4 )then
      
      call calcula_distancias_ao_hub(rede)
      
      open(999, file = 'fi_pQMF_N_10_a_7_RRN_deg_6_Plus_Star_998_alp'//trim(adjustl(alp_char2))//'.dat', status = 'unknown')
      do j1 = 1, rede%nodes
         write(999,*) rede%deg(j1), lista_distancias(j1), (f1(j1))**2.0d0
      enddo
 
 	  close(999)
      call system('tar -czf '//trim(adjustl('fi_pQMF_N_10_a_7_RRN_deg_6_Plus_Star_998_alp'//trim(adjustl(alp_char2))))//'.tar.gz '//trim(adjustl('fi_pQMF_N_10_a_7_RRN_deg_6_Plus_Star_998_alp'//trim(adjustl(alp_char2))))//'.dat')
      call system('rm '//trim(adjustl('fi_pQMF_N_10_a_7_RRN_deg_6_Plus_Star_998_alp'//trim(adjustl(alp_char2))))//'.dat')
      
      deallocate(f1)
      deallocate(lista_distancias)
   endif
      
   if(ind_lamb == 0) stop "Calculou o limiar. Nao vai rodar dinamica"
   
   lamb = lamb0 + 1.0d0 * ind_lamb * dlamb
     
   lambdaf = lamb
      
llbd: do while( lamb <= lambdaf )
      !============================================================
      ! tempo
      !============================================================
      write(lamb_char, '(f10.6)') lamb     
      n_it = 0
      n_sub = 0
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))//'_lamb_'//trim(adjustl(lamb_char))
         !==============================================================
         open(336, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_'))//'RRN_Plus_Star_tam_'//trim(adjustl(tam_char))//'_lamb_'//trim(adjustl(lamb_char))
         !==============================================================
         open(337, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
         !==============================================================
      if( .not. usouCopia)then
         
         rho0 = p0
         rho = rho0
         
         Y4_old = 1.0_dp / (1.0_dp * rede%nodes)         
         Y4 = Y4_old
         
         call condicao_inicial(rede, p0, 0.0d0)                  
         
         tempoEscrita = tempoEscrita0
      !=================================================================
      else
         usouCopia = .False.
      endif
      !=================================================================     
lt:   do while(t <= tf)
         
         !call tempoDeEscrever()
         
         write(336, *) t, rho
         write(337, *) t, Y4
         
         !##############################################################
         call rk4(dt, t, rede, rede%nodes, rede%sumDeg, alp, lamb, mu)
         !##############################################################
         
         rho = 0.0_dp
         
         Ii_neg = .False.
         Ii_gt_1 = .False.
         do i1 = 1, size(Ii)   ! Calcula rho
            if(lista_de_clusters(i1) /= i_comp_gigante) cycle
            if( Ii(i1) < 0.0_dp)then
               Ii_neg = .True.
               exit
            elseif( Ii(i1) > 1.0_dp )then
               Ii_gt_1 = .True.
               exit
            endif
            rho = rho + Ii(i1)
         enddo
         
         rho  = rho/( 1.0_dp * comp_gigante )
         
         sumIi2 = (sum(Ii**2.0_dp))**0.5_dp
         
		 Y4 = sum( (Ii/sumIi2)**4.0_dp )
         !##############################################################
         t = t + dt
         !##############################################################
         
         if( abs(rho - rho0)/rho0 < tole )then
            n_it = n_it + 1
            if( n_it >= int(2.0d0/dt) )then
               write(*,*) "Conv. lbd = ", lamb, " e rho = ", rho, "."
               exit lt
            endif
         elseif( abs(Y4 - Y4_old)/Y4_old < tole )then
            n_it = n_it + 1
            if( n_it >= int(2.0d0/dt) )then
               write(*,*) "Conv. lbd = ", lamb, " e Y4 = ", Y4, "."
               exit lt
            endif
            n_it = 0
         endif

!         if( ( abs(rho - rho0)/rho0 > 5.0d0 * tole ) .and. ( rho < 1.0_dp/(10.0_dp * rede%nodes) ) )then
!            n_sub = n_sub + 1
!            if ( n_sub >= int(5.0d0/dt) )then
!               write(*,*) "Abs subcr com lambda = ", lamb, " e rho = ", rho, "."              
!               exit lt
!            endif
!         else
!            n_sub = 0
!         endif
         
         !==============================================================
         !   Se deu ruim por causa de dt grande, a gente pega uma
         !   configuracao passada e abaixa o dt
         !   
         !==============================================================
         if( ( Ii_neg == .True. ) .or. ( Ii_gt_1 == .True. ) )then                 
            write(*,*) "Atualizou dt de dt = ", dt, " para dt = ", dt/2.0d0            
            usouCopia = .False.
            
            !call TestaSeTemCopiaEstadoSistema(usouCopia)
            
            dt = dt/2.0d0
            
            if( .not. usouCopia)then 
               rho0 = p0
         
               call condicao_inicial(rede, p0, 0.0d0)                  
         
               tempoEscrita = tempoEscrita0
               !=========================================================
            else   
               usouCopia = .False.
            endif            
            
            cycle lt
         endif   
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
      if( n_it >= int(2.0d0/dt) )then

         write(333, *) lamb, Y4
         write(334,*) lamb, rho   
         write(335,*) lamb, t
                     
      endif
      !=================================================================
      lamb = lamb + dlamb
      !================================================================= 
      call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia*')) )          
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

   !====================================================================
      subroutine leAtehOfimUnidimensional(label, nomeArquivo, var1, teveLeitura)      
         integer, intent(in) :: label
         character(len=*) :: nomeArquivo
         real(dp), intent(out) :: var1
         logical, intent(inout) :: teveLeitura
         integer :: st, res
         
         inquire( file=trim(adjustl(nomeArquivo)), exist=res ) 
         if(res)then
            open(label, file=trim(adjustl(nomeArquivo)), status='old')
            write(*,*) "Arquivo jah existe"
            do
               read(label, *, iostat = st) var1
               write(*,*) st
               if( st /= 0) exit
               teveLeitura = .True.
            enddo
         else
            open(label, file=trim(adjustl(nomeArquivo)), status='new')
            write(*,*) "Arquivo teve que ser criado"
         endif
   !====================================================================         
      end subroutine


   !====================================================================
      subroutine leAtehOfimBidimensional(label, nomeArquivo, var1, var2, teveLeitura)      
         integer, intent(in) :: label
         character(len=*) :: nomeArquivo
         real(dp), intent(out) :: var1, var2
         logical, intent(inout) :: teveLeitura
         integer :: st, res
         
         inquire( file=trim(adjustl(nomeArquivo)), exist=res ) 
         if(res)then
            open(label, file=trim(adjustl(nomeArquivo)), status='old')
            write(*,*) "Arquivo jah existe"
            do
               read(label, *, iostat = st) var1, var2               
               if( st /= 0) exit
               teveLeitura = .True.
            enddo
         else
            open(label, file=trim(adjustl(nomeArquivo)), status='new')
            write(*,*) "Arquivo teve que ser criado"
         endif
   !====================================================================         
      end subroutine

 
      subroutine criaRedeEClassificaClusters(this, tam_rede1, grau_RRN1, grau_max1)
        type(grafoRRN_Plus_Star) :: this
        integer :: tam_rede1
        integer :: grau_RRN1
        integer :: grau_max1
!#######################################################################   
        call this%ligaRRNStar(tam_rede1, grau_RRN1, grau_max1) 
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
            read(buffer,*) grau_RRN

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_max

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
            write(*,*) "Recebeu 10 arqumentos"
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
            read(buffer,*) grau_RRN

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_max

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
         dlamb = 0.0125_dp/(1.0_dp * divisor)
!#######################################################################
         p0 = 10.0_dp/(1.0_dp * tam_rede) !1.0d-2
!#######################################################################
!  Indice da semente
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
   subroutine kNN_e_clustering(this)
      class(grafo) :: this

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'tam_'//trim(adjustl(tam_char))//'.dat'
      call calcula_k_nn(this, .True., 800, trim(adjustl(arq_1)))
      close(800)

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering_'))//'tam_'//trim(adjustl(tam_char))//'.dat'
      call clustering(this,.True., 800, trim(adjustl(arq_1)))
      close(800)
   end subroutine
!=======================================================================   
    
end program
