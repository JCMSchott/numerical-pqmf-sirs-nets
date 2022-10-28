 module pqmf
   use types
   use geraRede
   use mod_tools
   implicit none
!#######################################################################   
   real(dp), allocatable :: I_i(:), R_i(:) ! S(:) = 1 - I(:) - R(:)
   real(dp), allocatable :: RI_ij(:), SI_ij(:), RR_ij(:)
!#######################################################################
   integer, allocatable :: st_spec(:)
!#######################################################################   
   character(len=500) :: local
   character(len=1000) :: local_arquivo, I_i_arquivo, S_i_arquivo
   character(len=1000) :: SI_ij_arquivo, RI_ij_arquivo, SS_ij_arquivo
   character(len=1000) :: lb_vs_I_i_arquivo, t_vs_I_im_arquivo
!#######################################################################
   real(dp) :: lambda
   real(dp) :: t
!#######################################################################   
   logical :: foi_I_i, foi_SI_ij, foi_RI_ij
   real(dp):: qualvalor
!#######################################################################   
   integer  :: nlambda
   real(dp) :: lambda0
   real(dp) :: lambdaf
   real(dp) :: Iest
!#######################################################################      
   real(dp) :: dlambda  
!#######################################################################  
   real(dp), allocatable :: k1_I_i(:), k2_I_i(:), k3_I_i(:), k4_I_i(:)
   real(dp), allocatable :: k1_R_i(:), k2_R_i(:), k3_R_i(:), k4_R_i(:)
!#######################################################################      
   real(dp), allocatable :: k1_RI_ij(:), k2_RI_ij(:), k3_RI_ij(:), k4_RI_ij(:)
   real(dp), allocatable :: k1_SI_ij(:), k2_SI_ij(:), k3_SI_ij(:), k4_SI_ij(:)
   real(dp), allocatable :: k1_RR_ij(:), k2_RR_ij(:), k3_RR_ij(:), k4_RR_ij(:)
!#######################################################################
   real(dp) :: argI_i(:), argR_i(:), argRI_ij(:), argSI_ij(:), argRR_ij(:)  
!#######################################################################
   contains

subroutine aloca_variaveis_dinamicas(this)
   class(grafo), intent(in) :: this
   integer :: sum_deg
   
   sum_deg = sum(this%deg)
   
      !######AlocaListas#####################
      if(allocated(I_i)) deallocate(I_i)
      allocate(I_i(this%nodes))
      !######################################
      if(allocated(R_i)) deallocate(R_i)
      allocate(R_i(this%nodes))
      !######################################
      if(allocated(RI_ij)) deallocate(RI_ij)
      allocate(RI_ij(sum_deg))
      !######################################
      if(allocated(SI_ij)) deallocate(SI_ij)
      allocate(SI_ij(sum_deg))
      !######################################
      if(allocated(RR_ij)) deallocate(RR_ij)
      allocate(RR_ij(sum_deg))
      !######################################
      if(allocated(st_spec)) deallocate(st_spec)
      allocate(st_spec(sum_deg))
      !######################################
      
      if(allocated(k1_I_i))then
         deallocate(k1_I_i) 
         deallocate(k2_I_i)
         deallocate(k3_I_i)
         deallocate(k4_I_i)
         allocate(k1_I_i(this%nodes))
         allocate(k2_I_i(this%nodes))
         allocate(k3_I_i(this%nodes))
         allocate(k4_I_i(this%nodes))
      endif

      if(allocated(k1_R_i))then
         deallocate(k1_R_i) 
         deallocate(k2_R_i)
         deallocate(k3_R_i)
         deallocate(k4_R_i)
         allocate(k1_R_i(this%nodes))
         allocate(k2_R_i(this%nodes))
         allocate(k3_R_i(this%nodes))
         allocate(k4_R_i(this%nodes))
      endif

      if(allocated(k1_RI_ij)then
         deallocate(k1_RI_ij)
         deallocate(k2_RI_ij)
         deallocate(k3_RI_ij)
         deallocate(k4_RI_ij)
         allocate(k1_RI_ij(sum_deg))
         allocate(k2_RI_ij(sum_deg))
         allocate(k3_RI_ij(sum_deg))
         allocate(k4_RI_ij(sum_deg))
      endif

      if(allocated(k1_SI_ij)then
         deallocate(k1_SI_ij)
         deallocate(k2_SI_ij)
         deallocate(k3_SI_ij)
         deallocate(k4_SI_ij)
         allocate(k1_SI_ij(sum_deg))
         allocate(k2_SI_ij(sum_deg))
         allocate(k3_SI_ij(sum_deg))
         allocate(k4_SI_ij(sum_deg))
      endif      

      if(allocated(k1_RR_ij)then
         deallocate(k1_RR_ij)
         deallocate(k2_RR_ij)
         deallocate(k3_RR_ij)
         deallocate(k4_RR_ij)
         allocate(k1_RR_ij(sum_deg))
         allocate(k2_RR_ij(sum_deg))
         allocate(k3_RR_ij(sum_deg))
         allocate(k4_RR_ij(sum_deg))
      endif            

!      real(dp) :: argI_i(this%nodes), argR_i(this%nodes), argRI_ij(sum_deg), argSI_ij(sum_deg), argSS_ij(sum_deg)
      
end subroutine      

   subroutine condicao_inicial(this)     

      class(grafo), intent(in) :: this
      integer:: j1, j12, j21, j23, j2, j3
      real(dp), parameter :: I_0 = 0.01_dp
      logical :: isI_i
      integer :: iost
      character(len=30) :: Nichar
      character(len=10) :: lambdachar
      character(len=100) :: arquivo_char2
      integer :: sum_deg
      
      sum_deg = sum(this%deg)
            
      foi_I_i = .False.
      foi_SI_ij = .False.
      foi_RI_ij = .False.
            
      !######Seta_Valor_Inicial#####################
      RI_ij = 0.0_dp
      !######AlocaListas############################
      t = 0.0_dp
                     
         do j2 = 1, this%nodes
            if(lista_de_clusters(j2) /= i_comp_gigante)then
               I_i(j2) = 0.0_dp
               R_i(j2) = 1.0_dp
               cycle
            endif
            I_i(j2) = I_0
            R_i(j2) = 0d0

            do j12 = this%aux(j2), this%aux(j2) + this%deg(j2) - 1
               !if(SI_ij(j12) > 0.0_dp) cycle
               j3 = this%listAdj(j12)
 
               !j21 = st_spec(j12)
                                        
               SI_ij(j12) = (1.d0 - I_i(j2)) * I_i(j3)
               
               !SI_ij(j21) = S_i(j3) * I_i(j2)
 
!               RR_ij(j12) = S_i(j2) * S_i(j3)
            
!               !SS_ij(j21) = SS_ij(j12)
                                 
            enddo            
         enddo                         
   end subroutine       
     
   !subroutine k4_sirs_pqmf(dt, t, this, ms, ns, alp, lbd, mu, I_i, S_i, RI_ij, SI_ij, SS_ij, fI_i, fS_i, fRI_ij, fSI_ij, fSS_ij)
   subroutine k4_sirs_pqmf(dt, t, this, ms, ns, alp, lbd, mu)      
      use types
      use geraRede
      
      real(dp) :: dt
      real(dp) :: t
      class(grafo), intent(in) :: this         
      integer :: ms
      integer :: ns
      real(dp) :: alp, lbd, mu

      integer :: l1
      
      !#########################SUBROTINA#############################

!#################################################################################################
		k1_I_i = dt * fI_i(t, this, ms, ns, lbd, mu, I_i, SI_ij)		
!#######################################################################		
		k1_S_i = dt * fS_i(t, this, ms, ns, alp, lbd, I_i, S_i, SI_ij)
!#######################################################################		
		k1_RI_ij = dt * fRI_ij(t, this, ms, ns, alp, lbd, mu, I_i, S_i, RI_ij, SI_ij, SS_ij)
!#######################################################################		
		k1_SI_ij = dt * fSI_ij(t, this, ms, ns, alp, lbd, mu, S_i, RI_ij, SI_ij, SS_ij)
!#######################################################################		
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
!#######################################################################						
		k2_S_i = dt * fS_i(t + (dt/2_dp), this, ms, ns, alp, lbd, argI_i, argS_i, argSI_ij)
!#######################################################################		
		k2_RI_ij = dt * fRI_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
!#######################################################################				
		k2_SI_ij = dt * fSI_ij(t + (dt/2_dp), this, ms, ns, alp, lbd, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
!#######################################################################				                
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
               I_i = I_i + 1d0/6d0 * (k1_I_i + 2d0 * k2_I_i + 2d0 * k3_I_i + k4_I_i)

               S_i = S_i + 1d0/6d0 * (k1_S_i + 2d0 * k2_S_i + 2d0 * k3_S_i + k4_S_i)

               RI_ij = RI_ij + 1d0/6d0 * (k1_RI_ij + 2d0 * k2_RI_ij + 2d0 * k3_RI_ij + k4_RI_ij)

               SI_ij = SI_ij + 1d0/6d0 * (k1_SI_ij + 2d0 * k2_SI_ij + 2d0 * k3_SI_ij + k4_SI_ij)
               
               SS_ij = SS_ij + 1d0/6d0 * (k1_SS_ij + 2d0 * k2_SS_ij + 2d0 * k3_SS_ij + k4_SS_ij)
!#########################SUBROTINA##############################################################
                                             
   end subroutine

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
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
!#######################################################################
         fI_i(k1) = -mu * I_i(k1)          
!#######################################################################          
         soma = 0.0_dp          
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            soma = soma + SI_ij(k12)
         enddo
!#######################################################################         
         fI_i(k1) = fI_i(k1) + lambda * soma
!#######################################################################
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
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
!#######################################################################         
         fS_i(k1) = alp * (1.0_dp -S_i(k1) -I_i(k1))
!#######################################################################
         soma = 0.0_dp          
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
            soma = soma + SI_ij(k12)
         enddo
!#######################################################################         
         fS_i(k1) = fS_i(k1) - lambda * soma    
!#######################################################################
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
      !#################################################################
      ! 	PARA TESTES
      real(dp) :: dum_II_ij2
      !#################################################################       
      integer :: edge_ji
      integer :: viz1
      integer :: k3    ! k1, k2 ja existem no geraRede
      integer :: k12, k12_dum, k21, k23, k32  
      real(dp) :: soma

      do k1 = 1, this%nodes
        if(lista_de_clusters(k1) /= i_comp_gigante) cycle
        ! Stubs que saem do sitio k1
        do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
           ! O vizinho do sitio k1
           k2 = this%listAdj(k12)           
!#######################################################################           
           k21 = st_spec(k12)
!#######################################################################
! Toda a fonte de erro esta nesse pedaco de codigo.
! Aparentemente, nao posso trocar i por j (1 por 2) nas duas linhas nao
! comentadas abaixo. Por que?
!#######################################################################
! Esse aqui tava dando errado. Mas to fazendo um teste
      dum_II_ij = I_i(k1) - SI_ij(k21) - RI_ij(k21)
!#######################################################################
! Esse aqui eh o certo. Depois do teste, reativo.           
      !dum_II_ij = I_i(k2) - SI_ij(k12) - RI_ij(k12)
!#######################################################################
      !dum_II_ij2 = I_i(k2) - SI_ij(k12) - RI_ij(k12)
!#######################################################################
!
!     if(abs(dum_II_ij-dum_II_ij2)>0.0_dp)then
!       write(*,*) "abs(dum_II_ij-dum_II_ij2)>0 e = ", abs(dum_II_ij-dum_II_ij2)
!       stop 
!     endif

 !    if(abs(SS_ij(k12)-SS_ij(k21))>0.0_dp)then
 !      write(*,*) "abs(SS_ij-SS_II_ji)>0 e = ", abs(SS_ij(k12)-SS_ij(k21))
 !      stop 
 !    endif

!      write(*,*) "II_ij = ", SS_ij(k12)
!      write(*,*) "II_ji = ", SS_ij(k21)

!      write(*,*) "II_ij = ", dum_II_ij                    !Eles sao diferentes, mas pq???
!      write(*,*) "II_ji = ", dum_II_ij2

!            write(*,*) "II_ij = ", I_i(k1) - SI_ij(k21) - RI_ij(k21)
!            write(*,*) "II_ji = ", I_i(k2) - SI_ij(k12) - RI_ij(k12)           !Eles dao absurdamente diferente aqui!
 

! Se II_ij eh simetrico, vale:

!       dum_II_ij = 0.5_dp * (I_i(k1) + I_i(k2) - SI_ij(k12) - RI_ij(k12) - SI_ij(k21) - RI_ij(k21))
!#######################################################################
! Mas, como o primeiro esta dando errado, esse tambem dara.
! O erro deve estar em SI ou RI, penso eu.
! I_i eh canonico, provavelmente nao eh nele.
           
!#######################################################################          

           fRI_ij(k12) = -(alp + mu) * RI_ij(k12) + mu * dum_II_ij
                      
           soma = 0.0_dp
           do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
              if(k23 == k21) cycle
              soma = soma + SI_ij(k23)
           enddo

! Esse aqui tava dando errado, mas to fazendo teste
! Nao! Nao era aqui o erro. Ele eh valido tambem.
            
           dum_RS_ijSj = (S_i(k2) -SI_ij(k21) - SS_ij(k12))/S_i(k2)   !anterior

           !Esse que tava dando certo. Depois do teste, volto a usar
           !dum_RS_ijSj = (S_i(k2) -SI_ij(k21) - SS_ij(k21))/S_i(k2)                      
 
           fRI_ij(k12) = fRI_ij(k12) + lambda * dum_RS_ijSj * soma
            
        enddo
      enddo
      
!     return
      
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
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle     
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) - 1 !a1
            
            k21 = st_spec(k12)
            
            k2 = this%listAdj(k12)            

            fSI_ij(k12) = -(mu + lambda) * SI_ij(k12) + alp * RI_ij(k12) 
            
            
            soma = 0.0_dp
            
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1               
               if(k23 == k21)cycle
               soma = soma + SI_ij(k23)
            enddo

            dum_SS_ijSj = SS_ij(k12)/S_i(k2)
            
            fSI_ij(k12) = fSI_ij(k12) + lambda * dum_SS_ijSj * soma
            
            !fSI_ij(k12) = fSI_ij(k12) + lambda * (SS_ij(k12)/S_i(k2)) * soma            
            
            
            soma = 0.0_dp
            
            do k12_dum = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
               if(k12_dum == k12) cycle
               soma = soma + SI_ij(k12_dum)
            enddo

             dum_SI_ijSi = SI_ij(k12)/S_i(k1)
             fSI_ij(k12) = fSI_ij(k12) - lambda * dum_SI_ijSi * soma
             
             !fSI_ij(k12) = fSI_ij(k12) - lambda * (SI_ij(k12)/S_i(k1)) * soma

         enddo
         
      enddo
      
!     return
      
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
         if(lista_de_clusters(k1) /= i_comp_gigante) cycle
         do k12 = this%aux(k1), this%aux(k1) + this%deg(k1) -1
            k2 = this%listAdj(k12)
       
            k21 = st_spec(k12)
                             
! Esse supostamente ta certo
           dum_RS_ij = S_i(k2) - SI_ij(k21) - SS_ij(k21)
           dum_RS_ji = S_i(k1) - SI_ij(k12) - SS_ij(k12)

! Esse supostamente nao ta certo
!			dum_RS_ij = S_i(k1) - SI_ij(k12) - SS_ij(k21)
!           dum_RS_ji = S_i(k2) - SI_ij(k21) - SS_ij(k12)
            
            fSS_ij(k12) = alp * (dum_RS_ij + dum_RS_ji)
            
            
            soma = 0.0_dp
            do k12_dum = this%aux(k1), this%aux(k1) + this%deg(k1) - 1
               if(k12_dum == k12) cycle
               soma = soma + SI_ij(k12_dum)
            enddo

            dum_SS_ijSi = SS_ij(k12)/S_i(k1)

            fSS_ij(k12) = fSS_ij(k12) - lambda * dum_SS_ijSi * soma
            

!            write(*,*) "II_ij = ", I_i(k1) - SI_ij(k21) - RI_ij(k21)
!            write(*,*) "II_ji = ", I_i(k2) - SI_ij(k12) - RI_ij(k12)           !Eles dao absurdamente diferente aqui!
            
            soma = 0.0_dp
            
            do k23 = this%aux(k2), this%aux(k2) + this%deg(k2) - 1
               if(k23 == k21) cycle
               soma = soma + SI_ij(k23)
            enddo

            dum_SS_ijSj = SS_ij(k12)/S_i(k2)

            fSS_ij(k12) = fSS_ij(k12) - lambda * dum_SS_ijSj * soma             
         enddo
      
      enddo      
   end function
   
end module

!#######################


program sirs_pqmp_est

   
!#####preambulo######   

!  use mod_rndgen
   use geraRede
   use mod_tools
   use types
   use pqmf
   use mod_tictoc
   !#####preambulo######

   implicit none
      
   integer, parameter :: t0 = 0.0_dp
   real(dp), parameter :: tf = 1000.0_dp
   real(dp), parameter :: dt = 0.01_dp
   integer, parameter :: npts = int((tf - t0)/dt)
   
   
   real(dp) :: alp
   real(dp), parameter :: mu = 1.0_dp
   integer :: semente
   
!  type(rndgen) :: gerador
   
   integer :: i1, i2, i3, i4, i5, i6, i7
   !integer :: j1   
 !######################################################################
 ! REDE
   type(grafo_PL_UCM) :: rede
   !integer :: tam_rede
   integer, allocatable :: tam_rede(:)
   real(dp) :: exp_gama
   integer :: k_i
   real(dp) :: k_f
   integer :: sum_deg
   integer :: sum_deg_comp_gig
 !######################################################################     
 ! VARIAVEIS DINAMICAS
   real(dp) :: I_im0
   real(dp) :: R_im0
   real(dp) :: RI_ijm0
   real(dp) :: S_im0
   real(dp) :: SI_ijm0
   real(dp) :: SS_ijm0
 !######################################################################        
   real(dp) :: I_im
   real(dp) :: R_im
   real(dp) :: RI_ijm
   real(dp) :: S_im
   real(dp) :: SI_ijm
   real(dp) :: SS_ijm
 !######################################################################         
   integer :: I_im_repetido
   real(dp) :: delta_I_im
   logical :: tem_trans
 !######################################################################     
   character(len=100) :: arquivo_char1, arquivo_rede
   real(dp), parameter :: t_rec = 1.0_dp * int(tf/10.0_dp)
   real(dp) :: t_rec1
   logical :: arquivo_existe
   integer :: n_args
   integer :: flag1
   character(len=1) :: flag_char
   character(len=10) :: lambda_char
   character(len=10) :: alp_char
   character(len=7) :: tam_char   
   character(len=7) :: gama_char
   character(len=300) :: buffer
   integer :: label
   real(dp) :: erro
   real(dp), parameter :: toll = 1d-9
 !######################################################################
   integer(kind=8) :: t1, t2, taxa
   real(dp) :: te 
   integer :: m1
 !######################################################################
   type(tictoc) :: crono
 !######################################################################        
   local = trim(adjustl('/home/jota/SIRS_pQMF/PQMF_Rst_Rede_Sint/'))   
 !#######################################################################
!  ESCOLHE SE VAI VOLTAR DE UMA CONFIGURACAO ANTIGA OU NAO
!#######################################################################
!   if(n_args == 1)then
!      call getarg(1, flag_char)
 !     read(flag_char,*) flag1
!   else
 !     stop "Numero de argumentos invalido"
 !  endif   
!####################################################################### 
! DEFINO PARAMETROS LAMBDA.
!  OBS.:

!#######################   
   nlambda = 1000
   dlambda = 0.00125_dp   
   lambda0 = dlambda + 0.09_dp ! +0.5_dp
   lambda = lambda0
   lambdaf = lambda0 + 1.0_dp * nlambda * dlambda
!#######################   
   alp = 1.5_dp  
!################################################################################   
! DEFINO SEMENTE, GRAU MINIMO, TAMANHO DA REDE E EXP_GAMA
!#######################
   semente  = 967891968   
!#######################
   k_i = 3   
!#######################   
   !tam_rede = 1000
!#######################   
   exp_gama = 2.3_dp
   !exp_gama = 2.7_dp
   !exp_gama = 3.5_dp
!##############################################################################################
allocate(tam_rede(9))

tam_rede(1) = 1000; tam_rede(2) = 3000

tam_rede(3) = 10000; tam_rede(4) = 30000

tam_rede(5) = 100000; tam_rede(6) = 300000

tam_rede(7) = 1000000; tam_rede(8) = 3000000
         
!#############################################################################################
!#######################
   write(gama_char, '(f5.2)') exp_gama         
   gama_char = trim(adjustl(gama_char))
!#######################
!##########################################################################################

arquivo_char1 = trim(adjustl('N_vs_te_rk4'//'_gam_'//trim(adjustl(gama_char))))

local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'


open(unit = 888, file = trim(adjustl(local_arquivo)), status = 'unknown')
!##########################################################################################

call crono%start()

tmn:  do m1 = 1, size(tam_rede)

  
! ESCOLHE A STRING CORRESPONDENTE AO TAMANHO DA REDE  
!if(.True.)then 
!   if(tam_rede == 100)then
!      tam_char = '0_1k'
!   elseif(tam_rede == 1000)then
!      tam_char = '1k'
!   elseif(tam_rede == 3000)then
!      tam_char = '3k'
!   elseif(tam_rede == 5000)then
!      tam_char = '5k'      
!   elseif(tam_rede == 10000)then
!      tam_char = '10k'
!   elseif(tam_rede == 30000)then
!      tam_char = '30k'
!   elseif(tam_rede == 100000)then
!      tam_char = '100k'
!   elseif(tam_rede == 300000)then
!      tam_char = '300k'
!   elseif(tam_rede == 1000000)then
!      tam_char = '1M'
!   elseif(tam_rede == 3000000)then
!      tam_char = '3M'
!   elseif(tam_rede == 10000000)then
!      tam_char = '10M'
!   elseif(tam_rede == 30000000)then
!      tam_char = '30M'
!   else
!      stop 'Escolha um tamanho de rede dentro do catalogo'
!   endif
!endif   


! ESCOLHE A STRING CORRESPONDENTE AO TAMANHO DA REDE  
if(.True.)then 
   if(tam_rede(m1) == 100)then
      tam_char = '0_1k'
   elseif(tam_rede(m1) == 1000)then
      tam_char = '1k'
   elseif(tam_rede(m1) == 3000)then
      tam_char = '3k'
   elseif(tam_rede(m1) == 5000)then
      tam_char = '5k'      
   elseif(tam_rede(m1) == 10000)then
      tam_char = '10k'
   elseif(tam_rede(m1) == 30000)then
      tam_char = '30k'
   elseif(tam_rede(m1) == 100000)then
      tam_char = '100k'
   elseif(tam_rede(m1) == 300000)then
      tam_char = '300k'
   elseif(tam_rede(m1) == 1000000)then
      tam_char = '1M'
   elseif(tam_rede(m1) == 3000000)then
      tam_char = '3M'
   elseif(tam_rede(m1) == 10000000)then
      tam_char = '10M'
   elseif(tam_rede(m1) == 30000000)then
      tam_char = '30M'
   else
      stop 'Escolha um tamanho de rede dentro do catalogo'
   endif
endif   



!##########################################################################################
!				Inicia grafo
!##########################################################################################
if(.True.)then
!   if(exp_gama > 3.0_dp)then
!      k_f = 1.0_dp * (1.0_dp * tam_rede)**(1.0_dp/(exp_gama -1.0_dp))
!   else
!      k_f = 1.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
!   endif      

   if(exp_gama > 3.0_dp)then
      k_f = 1.0_dp * (1.0_dp * tam_rede(m1))**(1.0_dp/(exp_gama -1.0_dp))
   else
      k_f = 1.0_dp * (1.0_dp * tam_rede(m1))**(0.5_dp)
   endif      
  
!##########################################################################################         
   !call rede%iniciaGrafo(tam_rede)
   call rede%iniciaGrafo(tam_rede(m1))
   call rede%inicia(k_i, k_f, exp_gama, semente)   
   call rede%liga(semente, .True.)
            
   sum_deg = sum(rede%deg)
endif   
!##########################################################################################         
! DEFINE AS STRINGS
!#######################           
   tam_char = trim(adjustl(tam_char))
!#######################
   write(alp_char, '(f5.2)') alp
   alp_char = trim(adjustl(alp_char))   
!#######################################################################   
   call sub_classifica_clusters(rede, .False., 000, 'nenhum.dat')
!#######################################################################
   !write(tam_char,*) comp_gigante
   
   sum_deg = sum(rede%deg)
        
   sum_deg_comp_gig = 0
   
   call aloca_variaveis_dinamicas(rede)
            
   do i1 = 1, rede%nodes 
      if(lista_de_clusters(i1) /= i_comp_gigante) cycle      
      do i2 = rede%aux(i1), rede%aux(i1) + rede%deg(i1) -1
         do i3 = rede%aux(rede%listAdj(i2)), rede%aux(rede%listAdj(i2)) + rede%deg(rede%listAdj(i2)) - 1
            if(rede%listAdj(i3) /= i1) cycle
            st_spec(i2) = i3
         enddo 
      enddo      
      sum_deg_comp_gig = sum_deg_comp_gig + rede%deg(i1)
   enddo
            
   write(*,*) "Calculou componente gigante"
   write(*,*) ""    

 write(*,*) "##########################################################"  
 write(*,*) "Instante de tempo inicial t0 = ",t0
 write(*,*) "Instante de tempo final tf = ", tf
 write(*,*) "Incremento de tempo dt = ", dt
 write(*,*) "Numero de passos de tempo = ", npts
 write(*,*) "##########################################################"   
 write(*,*) "Taxa de infeccao inicial lambda0 = ", lambda0
 write(*,*) "Incremento da taxa de infeccao dlambda = ", dlambda
 write(*,*) "Numero de passos de taxa de infeccao = ", nlambda
 write(*,*) "Taxa de enfraquecimento imunologico alfa = ", alp   
 write(*,*) "##########################################################"
 write(*,*) "Tamanho da rede = ", rede%nodes
 write(*,*) "Expoente gama = ", exp_gama
 write(*,*) "##########################################################"    
 
   
if(.True.)then ! ABRIA UM MONTE DE ARQUIVO PRA CHECAR DINAMICA           
!############################################################################################
!   arquivo_char1 = trim(adjustl('t_vs_I_i_tam_'//trim(adjustl(tam_char))//'_rede_'//trim(adjustl(arquivo_rede))))
         
!   local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
                  
!   I_i_arquivo = local_arquivo
         
!   open(unit = 26, file = trim(adjustl(I_i_arquivo)), status ='unknown')
                     
!############################################################################################
!   arquivo_char1 = trim(adjustl('t_vs_S_i_tam_'//trim(adjustl(tam_char))//'_rede_'//trim(adjustl(arquivo_rede))))

!   local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            
!         S_i_arquivo = local_arquivo
!         open(unit = 27, file = trim(adjustl(S_i_arquivo)), status ='unknown')
!############################################################################################
            
!         arquivo_char1 = trim(adjustl('t_vs_SI_ij_tam_'//trim(adjustl(tam_char))//'_rede_'//trim(adjustl(arquivo_rede))))
         
!         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            
!         SI_ij_arquivo = local_arquivo
            
!         open(unit = 28, file = trim(adjustl(SI_ij_arquivo)), status = 'unknown')
         
!############################################################################################
!         arquivo_char1 = trim(adjustl('t_vs_RI_ij_tam_'//trim(adjustl(tam_char))//'_rede_'//trim(adjustl(arquivo_rede))))
!         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
!         RI_ij_arquivo = local_arquivo
            
!         open(unit = 29, file = trim(adjustl(RI_ij_arquivo)), status = 'unknown')
         
!############################################################################################
!         arquivo_char1 = trim(adjustl('t_vs_SS_ij_tam_'//trim(adjustl(tam_char))//'_rede_'//trim(adjustl(arquivo_rede))))
            
!         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            
!         SS_ij_arquivo = local_arquivo
            
!         open(unit = 30, file = trim(adjustl(SS_ij_arquivo)), status = 'unknown')

!        endif
endif
 
 if(.True.)then !AQUI ABRO OS ARQUIVOS ESTACIONARIOS  
!############################################################################################
         arquivo_char1 = trim(adjustl('lbd_vs_I_im_tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))

         local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
         open(unit = 20, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
!############################################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_S_im_tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))

            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
            open(unit = 21, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
!############################################################################################

!############################################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_R_im_tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
         
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
            open(unit = 22, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')

!############################################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_SI_ijm_tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
         
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat' 
         
            open(unit = 23, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')

!############################################################################################

            arquivo_char1 = trim(adjustl('lbd_vs_RI_ijm_tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
         
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'

            open(unit = 24, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
!############################################################################################

            arquivo_char1 = trim(adjustl('lbd_vs_SS_ijm_tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))

            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
         
            open(unit = 25, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')

!############################################################################################

            arquivo_char1 = trim(adjustl('lbd_vs_t_conv_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))

            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'


            open(unit = 32, file = trim(adjustl(local_arquivo)), access = 'append', status = 'unknown')
endif                  
!#############################################################################################        
         do i3 = 1, nlambda
!#############################################################################################            
            call condicao_inicial(rede) 
!#############################################################################################                                    
            if(.True.)then ! ABRO ARQUIVOS PRA CHECAR CONVERGENCIA   
               write(lambda_char, '(f10.7)') lambda
               arquivo_char1 = trim(adjustl('t_vs_I_im_tam_'//trim(adjustl(tam_char))//'_alp_'//trim(adjustl(alp_char))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 31, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
!############################################################################################

               arquivo_char1 = trim(adjustl('t_vs_dI_im_I_im_tam_'//trim(adjustl(tam_char))//'_alp_'//trim(adjustl(alp_char))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 33, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
            endif               
!#######################################################################
! 	VALORES INICIAIS DAS VARIAVEIS DINAMICAS
!#######################################################################
			if(.True.)then ! SETTO VALORES INICIAIS DAS MEDIAS
				I_im = 0.0_dp
				do i4 = 1, size(I_i)
				   if(lista_de_clusters(i4) /= i_comp_gigante) cycle
				   I_im = I_im + I_i(i4)
				enddo
				I_im = 1.0_dp * I_im/comp_gigante
	!#######################################################################            
				S_im = 0.0_dp
				do i4 = 1, size(S_i)
				   if(lista_de_clusters(i4) /= i_comp_gigante) cycle
				   S_im = S_im + S_i(i4)
				enddo            
				S_im = 1.0_dp * S_im/comp_gigante !size(S_i)           
	!#######################################################################            
				R_im = 0.0_dp
				do i4 = 1, size(I_i)
				   if(lista_de_clusters(i4) /= i_comp_gigante) cycle
				   R_im = R_im + 1.0_dp - I_i(i4) - S_i(i4)
				enddo                        
				R_im = 1.0_dp - I_im - S_im            
	!#######################################################################            
				SI_ijm = 0.0_dp
				do i4 = 1, size(SI_ij)
				   if(lista_de_clusters(rede%listAdj(i4)) /= i_comp_gigante) cycle
				   SI_ijm = SI_ijm + SI_ij(i4)
				enddo                                    
				SI_ijm = 1.0_dp * SI_ijm/sum_deg_comp_gig !size(SI_ij)
	!#######################################################################
				SS_ijm = 0.0_dp
				do i4 = 1, size(SS_ij)
				   if(lista_de_clusters(rede%listAdj(i4)) /= i_comp_gigante) cycle
				   SS_ijm = SS_ijm + SS_ij(i4)
				enddo                        
				SS_ijm = 1.0_dp * SS_ijm/sum_deg_comp_gig !size(SS_ij)
	!#######################################################################
				RI_ijm = 0.0_dp
				do i4 = 1, size(RI_ij)
				   if(lista_de_clusters(rede%listAdj(i4)) /= i_comp_gigante) cycle
				   RI_ijm = RI_ijm + RI_ij(i4)
				enddo                        
				RI_ijm = 1.0_dp * RI_ijm/sum_deg_comp_gig !size(RI_ij)
			endif	            
	!#######################################################################
            I_im_repetido = 0
            
            t_rec1 = t_rec
                     
din1:       do i4 = 1, npts
                           
               I_im0 = I_im                              
               S_im0 = S_im
               R_im0 = 1.0_dp - S_im0 - I_im0
               
               RI_ijm0 = RI_ijm
               SI_ijm0 = SI_ijm
               SS_ijm0 = SS_ijm
               
               
               write(31,*) t, I_im
                   
               !call crono%reset()
               
               !call crono%tic()               
               
               !call k4_sirs_pqmf(dt, t, rede, rede%nodes, sum_deg, alp, lambda, mu, I_i, S_i, RI_ij, SI_ij, SS_ij, fI_i, fS_i, fRI_ij, fSI_ij, fSS_ij)
               
               !call crono%toc()
               
               call k4_sirs_pqmf(dt, t, rede, rede%nodes, sum_deg, alp, lambda, mu)
!#######################################################################
               
               !te = 1.0_dp *(t2-t1)/(1.0_dp * taxa)
               
               !write(888,*) rede%nodes, crono%t_tot
               
               !write(*,*) rede%nodes, crono%t_tot
               
               !stop
               cycle tmn
               
               I_im = 0.0_dp
               do i5 = 1, size(I_i)
                  if(lista_de_clusters(i5) /= i_comp_gigante) cycle
                  if(I_i(i5) < 0.0_dp)then
                     write(*,*) "I_i = ", I_i(i5)
                     stop
                  endif
                  I_im = I_im + I_i(i5)
               enddo
               I_im = 1.0_dp * I_im/comp_gigante
!#######################################################################
               S_im = 0.0_dp
               do i5 = 1, size(S_i)
                  if(lista_de_clusters(i5) /= i_comp_gigante) cycle
                  if(S_i(i5) < 0.0_dp)then
                     write(*,*) "S_i = ", S_i(i5)
                     stop
                  endif
                  S_im = S_im + S_i(i5)
               enddo
               S_im = 1.0_dp * S_im/comp_gigante
!#######################################################################
               R_im = 0.0_dp
               do i5 = 1, size(I_i)
                  if(lista_de_clusters(i5) /= i_comp_gigante) cycle

!#######################################################################
! FACO TESTES PRA SABER SE II_ij E II_ji SAO SIMETRICOS
!                  do i6 = rede%aux(i5), rede%aux(i5) + rede%deg(i5) - 1
!                     i7 = rede%listAdj(i6)
                     
!                     write(*,*) "II_ij, II_ji ", I_i(i5) - SI_ij(st_spec(i6)) - RI_ij(st_spec(i6)), I_i(i7) - SI_ij(i6) - RI_ij(i6)
!                     write(*,*) "II_ij, II_ji ", abs(I_i(i5) - SI_ij(st_spec(i6)) - RI_ij(st_spec(i6)) -(I_i(i7) - SI_ij(i6) - RI_ij(i6)))/(I_i(i7) - SI_ij(i6) - RI_ij(i6))                  
!                  enddo
!#######################################################################                  
                                    
                  if((1.0_dp - S_i(i5) - I_i(i5)) < 0.0_dp)then
                     write(*,*) "R_i = ", 1.0_dp - S_i(i5) - I_i(i5)
                     stop
                  endif
                  R_im = R_im + 1.0_dp - I_i(i5) - S_i(i5)
               enddo
               R_im = 1.0_dp * R_im/comp_gigante
!#######################################################################
            do i5 = 1, size(RI_ij)
               if(lista_de_clusters(rede%listAdj(i5)) /= i_comp_gigante) cycle
               if(RI_ij(i5) < 0.0_dp)then
                  write(*,*) "RI_ij = ", RI_ij(i5)
                  stop
               endif
               RI_ijm = RI_ijm + RI_ij(i5)
            enddo                        
            RI_ijm = 1.0_dp * RI_ijm/sum_deg_comp_gig
!#######################################################################
            do i5 = 1, size(SI_ij)
               if(lista_de_clusters(rede%listAdj(i5)) /= i_comp_gigante) cycle
               if(SI_ij(i5) < 0.0_dp)then
                  write(*,*) "SI_ij = ", SI_ij(i5)
                  stop
               endif
               SI_ijm = SI_ijm + SI_ij(i5)
            enddo                        
            SI_ijm = 1.0_dp * SI_ijm/sum_deg_comp_gig
!#######################################################################
            do i5 = 1, size(SS_ij)
               if(lista_de_clusters(rede%listAdj(i5)) /= i_comp_gigante) cycle
               
               !write(*,*) "SS_ij, SS_ji = ", SS_ij(i5), SS_ij(st_spec(i5))
               
               if(SS_ij(i5) > 1.0_dp)then
                  write(*,*) "SS_ij = ", SS_ij(i5)
                  stop
               endif
               SS_ijm = SS_ijm + SS_ij(i5)
            enddo                        
            SS_ijm = 1.0_dp * SS_ijm/sum_deg_comp_gig
!#######################################################################               
            write(33,*) t, abs(I_im - I_im0)/I_im0 !abs(I_im - I_im0)/I_im0
                                              
               if((abs(I_im - I_im0)/I_im0) < toll)then
                  I_im_repetido = I_im_repetido + 1
                  if(I_im_repetido == 10)then
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
               
               if(I_im < 1.0_dp/(rede%nodes))then
                 I_im = 0.0_dp
                  S_im = 1.0_dp
                  R_im = 0.0_dp
                  RI_ijm = 0.0_dp
                  SI_ijm = 0.0_dp
                  SS_ijm = 1.0_dp
                  exit din1
               endif
                              
               t = t + dt
                
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
enddo tmn
    
         !############Dinamica#################
         do i5 = 20, 32
            !if(i5 == 31)cycle
            close(i5)
         enddo
         close(33)   
!  stop
         close(888)         
         write(*,*) "Finalizou rotina"


end program
