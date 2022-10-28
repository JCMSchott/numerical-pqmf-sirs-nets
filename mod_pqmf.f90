module pqmf
   use types
   use geraRede
   implicit none
   
   real(dp), allocatable :: I_i(:), S_i(:) ! R(:) = 1 - I(:) - S(:)
   real(dp), allocatable :: RI_ij(:), SI_ij(:), SS_ij(:)
   
   contains
      
   subroutine condicao_inicial(this, semente)
      class(grafo), intent(in) :: this
      integer:: j1, j12, j21, j2, j23, j3
      integer :: semente
      integer :: sum_deg
      real(dp), allocatable :: II_ij(:)
      real(dp), allocatable :: RR_ij(:)
      real(dp), allocatable :: RS_ij(:)
      logical, allocatable :: Iaux_i(:)
      logical, allocatable :: Saux_i(:)
      logical, allocatable :: IIaux_ij(:)
      logical, allocatable :: SIaux_ij(:)
      logical, allocatable :: RIaux_ij(:)
      logical, allocatable :: RRaux_ij(:)
      logical, allocatable :: SSaux_ij(:)
      real(dp) :: sorteio
      type(rndgen) :: gen
      
      call gen%init(semente)
      
      sum_deg = sum(this%deg)
      
      !######AlocaListas#####################
      if(allocated(I_i)) deallocate(I_i)
      allocate(I_i(this%nodes))
      I_i = 0.0_dp

      if(allocated(Iaux_i)) deallocate(Iaux_i)
      allocate(Iaux_i(this%nodes))
      Iaux_i = .False.

      if(allocated(S_i)) deallocate(S_i)
      allocate(S_i(this%nodes))
      S_i = 0.0_dp

      if(allocated(II_ij)) deallocate(II_ij)
      allocate(II_ij(sum_deg))      
      II_ij = 0.0_dp

      if(allocated(IIaux_ij)) deallocate(IIaux_ij)
      allocate(IIaux_ij(sum_deg))      
      IIaux_ij = .False.

      if(allocated(RI_ij)) deallocate(RI_ij)
      allocate(RI_ij(sum_deg))      
      RI_ij = 0.0_dp

      if(allocated(RR_ij)) deallocate(RR_ij)
      allocate(RR_ij(sum_deg))      
      RR_ij = 0.0_dp

      if(allocated(RRaux_ij)) deallocate(RRaux_ij)
      allocate(RRaux_ij(sum_deg))      
      RRaux_ij = .False.

      if(allocated(RS_ij)) deallocate(RS_ij)
      allocate(RS_ij(sum_deg))      
      RS_ij = 0.0_dp

      if(allocated(RIaux_ij)) deallocate(RIaux_ij)
      allocate(RIaux_ij(sum_deg))      
      RIaux_ij = .False.

      if(allocated(SI_ij)) deallocate(SI_ij)
      allocate(SI_ij(sum_deg))      
      SI_ij = 0.0_dp

      if(allocated(SIaux_ij)) deallocate(SIaux_ij)
      allocate(SIaux_ij(sum_deg))      
      SIaux_ij = .False.

      if(allocated(SS_ij)) deallocate(SS_ij)
      allocate(SS_ij(sum_deg))
      SS_ij = 0.0_dp

      if(allocated(SSaux_ij)) deallocate(SSaux_ij)
      allocate(SSaux_ij(sum_deg))
      SSaux_ij = .False.
      
      !######AlocaListas#####################
      
      do j1 = 1, this%nodes
         
         if(.not. Iaux_i(j1))then
            I_i(j1) = gen%rnd()/10.0_dp
            Iaux_i(j1) = .True.
            
            S_i(j1) = (1.0_dp - I_i(j1)) * gen%rnd()/100.0_dp               
            
         endif
         
         do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) -1
            j2 = this%listAdj(j12)
                              
            do j23 = this%aux(j2), this%aux(j2) + this%deg(j2) -1
               if(this%listAdj(j23) == j1)then
                  j21 = j23
                  exit
               endif
            enddo
                              
            if(.not. IIaux_ij(j12))then
               II_ij(j12) = I_i(j1) * gen%rnd()/1000.0_dp 
               II_ij(j21) = II_ij(j12)
                  
               IIaux_ij(j21) = .True.
               IIaux_ij(j12) = .True.
                                              
            endif   

            if(.not. SIaux_ij(j21))then
                  
               SI_ij(j21) = (I_i(j1) - II_ij(j12)) * gen%rnd()/100.0_dp
               
!              if(SI_ij(j21) < 0.0_dp) stop "SI_ij(j21) < 0.0"
               
               RI_ij(j21) = I_i(j1) -(II_ij(j12) + SI_ij(j21))
               
!              if(RI_ij(j21) < 0.0_dp) stop "RI_ij(j21) < 0.0"
                  
               SIaux_ij(j21) = .True.
                  
            endif


            if(.not. SSaux_ij(j12))then
               SS_ij(j12) = S_i(j1) * gen%rnd()/1000.0_dp
               SS_ij(j21) = SS_ij(j12)
               
               SSaux_ij(j21) = .True.
               SSaux_ij(j12) = .True.
            
            endif            

            SI_ij(j12) = (S_i(j1) - SS_ij(j12)) * gen%rnd()/10.0_dp
            
            RS_ij(j21) = S_i(j1) -SI_ij(j12) - SS_ij(j12)
            
!           if(RS_ij(j21) < 0.0_dp) stop "RS_ij(j21) < 0.0"
            
            if(.not. RRaux_ij(j12))then
               RR_ij(j12) = (1.0_dp - I_i(j1) - S_i(j1)) * gen%rnd()/1000.0_dp
               RR_ij(j21) = RR_ij(j12)
               RRaux_ij(j12) = .True.
               RRaux_ij(j21) = .True.
   
            endif   
               
            RI_ij(j12) = (1.0_dp - I_i(j1) - S_i(j1) - RR_ij(j12)) * gen%rnd()/10.0_dp
            
!            if(RI_ij(j12) < 0.0_dp)then
!               RR_ij(j12) = (1.0_dp - I_i(j1) - S_i(j1)) * gen%rnd()/1000.0_dp
!               RR_ij(j21) = RR_ij(j12)
!               RI_ij(j12) = (1.0_dp - I_i(j1) - S_i(j1) - RR_ij(j12)) * gen%rnd()/10.0_dp
!            endif
            if(RI_ij(j12) < 0.0_dp) stop "RI_ij(j12) < 0.0"
            
            RS_ij(j12) = 1.0_dp - I_i(j1) - S_i(j1) - RI_ij(j12) - RR_ij(j12)
            
!           if(RS_ij(j12) < 0.0_dp) stop "RS_ij(j12) < 0.0"
            
            if(.not.Iaux_i(j2))then   
               I_i(j2) = II_ij(j21) + RI_ij(j12) + SI_ij(j12)
               S_i(j2) = SI_ij(j21) + SS_ij(j12) + RS_ij(j12)
               Iaux_i(j2) = .True.               
            endif                                    
         enddo
                      
      enddo

      !######DesalocaListas#####################      
      deallocate(Iaux_i)
            
      deallocate(II_ij)
      
      deallocate(IIaux_ij)
            
      deallocate(RR_ij)
      
      deallocate(RRaux_ij)
      
      deallocate(RS_ij)
      
      deallocate(RIaux_ij)
            
      deallocate(SIaux_ij)

      deallocate(SSaux_ij)
                   
   end subroutine       

     
   subroutine k4_sirs_pqmf(dt, t, this, m_sitios, n_arestas, alp, lambda, mu, I_i, S_i, RI_ij, SI_ij, SS_ij, fI_i, fS_i, fRI_ij, fSI_ij, fSS_ij)
      use types
      use geraRede
      
      real(dp) :: dt
      real(dp) :: t
      class(grafo), intent(in) :: this         
      integer :: m_sitios
      integer :: n_arestas
      real(dp) :: alp, lambda, mu

      real(dp) :: I_i(m_sitios)
      real(dp) :: S_i(m_sitios)
      real(dp) :: RI_ij(n_arestas)
      real(dp) :: SI_ij(n_arestas)
      real(dp) :: SS_ij(n_arestas)


      real(dp) :: I_ip(m_sitios)
      real(dp) :: S_ip(m_sitios)
      
      real(dp) :: k1_I_i(m_sitios), k2_I_i(m_sitios), k3_I_i(m_sitios), k4_I_i(m_sitios)
      real(dp) :: k1_S_i(m_sitios), k2_S_i(m_sitios), k3_S_i(m_sitios), k4_S_i(m_sitios)


      real(dp) :: RI_ijp(n_arestas)
      real(dp) :: SI_ijp(n_arestas)
      real(dp) :: SS_ijp(n_arestas)
      
      real(dp) :: k1_RI_ij(n_arestas), k2_RI_ij(n_arestas), k3_RI_ij(n_arestas), k4_RI_ij(n_arestas)
      real(dp) :: k1_SI_ij(n_arestas), k2_SI_ij(n_arestas), k3_SI_ij(n_arestas), k4_SI_ij(n_arestas)
      real(dp) :: k1_SS_ij(n_arestas), k2_SS_ij(n_arestas), k3_SS_ij(n_arestas), k4_SS_ij(n_arestas)

      real(dp) :: argI_i(m_sitios), argS_i(m_sitios), argRI_ij(n_arestas), argSI_ij(n_arestas), argSS_ij(n_arestas)


      !#########################INTERFACE#############################
      
      interface
 	function fI_i(t, this, m_sitios, n_arestas, lambda, mu, I_i, SI_ij)
           use types
           use geraRede
            
            real(dp) :: t
            class(grafo), intent(in) :: this       
            integer :: m_sitios
            integer :: n_arestas
            real(dp) :: lambda, mu
            real(dp) :: I_i(m_sitios)
            real(dp) :: SI_ij(n_arestas)       
            real(dp) :: fI_i(m_sitios) 
 	end function      


 	function fS_i(t, this, m_sitios, n_arestas, alp, lambda, I_i, S_i, SI_ij)
            use types
            use geraRede
            
            real(dp) :: t
            class(grafo), intent(in) :: this         
            integer :: m_sitios
            integer :: n_arestas
            real(dp) :: alp, lambda
            real(dp) :: I_i(m_sitios)
            real(dp) :: S_i(m_sitios)
            real(dp) :: SI_ij(n_arestas)       
            real(dp) :: fS_i(m_sitios)
 	end function 			
 			
        function fRI_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, I_i, S_i, RI_ij, SI_ij, SS_ij)
            use types
            use geraRede

            real(dp) :: t
            class(grafo), intent(in) :: this         
            integer :: m_sitios
            integer :: n_arestas
            real(dp) :: alp, lambda, mu
            real(dp) :: I_i(m_sitios)
            real(dp) :: S_i(m_sitios)
            real(dp) :: RI_ij(n_arestas)
            real(dp) :: SI_ij(n_arestas)
            real(dp) :: SS_ij(n_arestas)           
            real(dp) :: fRI_ij(n_arestas)
         end function  

         function fSI_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, S_i, RI_ij, SI_ij, SS_ij)

            use types
            use geraRede

            real(dp) :: t
            class(grafo), intent(in) :: this         
            integer :: m_sitios
            integer :: n_arestas
            real(dp) :: alp, lambda, mu
            real(dp) :: S_i(m_sitios)
            real(dp) :: RI_ij(n_arestas)
            real(dp) :: SI_ij(n_arestas)
            real(dp) :: SS_ij(n_arestas)           
            real(dp) :: fSI_ij(n_arestas) 			
 			end function

         function fSS_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, S_i, SI_ij, SS_ij)

            use types
            use geraRede

            real(dp) :: t
            class(grafo), intent(in) :: this         
            integer :: m_sitios
            integer :: n_arestas
            real(dp) :: alp, lambda, mu
            real(dp) :: S_i(m_sitios)
            real(dp) :: SI_ij(n_arestas)
            real(dp) :: SS_ij(n_arestas)           
            real(dp) :: fSS_ij(n_arestas)
 			end function 			      
      end interface        
      !#########################INTERFACE#############################


      
      !#########################SUBROTINA#############################

  	!####################################################################################
				
		k1_I_i = dt * fI_i(t, this, m_sitios, n_arestas, lambda, mu, I_i, SI_ij)		
		
		k1_S_i = dt * fS_i(t, this, m_sitios, n_arestas, alp, lambda, I_i, S_i, SI_ij)
		
		k1_RI_ij = dt * fRI_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, I_i, S_i, RI_ij, SI_ij, SS_ij)
		
		k1_SI_ij = dt * fSI_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, S_i, RI_ij, SI_ij, SS_ij)
		
		k1_SS_ij = dt * fSS_ij(t, this, m_sitios, n_arestas, alp, lambda, mu, S_i, SI_ij, SS_ij)
		
		!###############################		
		argI_i = I_i + 1d0/2d0 * k1_I_i
		argS_i = S_i + 1d0/2d0 * k1_S_i
		argRI_ij = RI_ij + 1d0/2d0 * k1_RI_ij
                argSI_ij = SI_ij + 1d0/2d0 * k1_SI_ij			
		argSS_ij = SS_ij + 1d0/2d0 * k1_SS_ij
		!###############################

		k2_I_i = dt * fI_i(t + (dt/2d0), this, m_sitios, n_arestas, lambda, mu, argI_i, argSI_ij)		
		
		k2_S_i = dt * fS_i(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, argI_i, argS_i, argSI_ij)
		
		k2_RI_ij = dt * fRI_ij(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k2_SI_ij = dt * fSI_ij(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k2_SS_ij = dt * fSS_ij(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, mu, argS_i, argSI_ij, argSS_ij)


		!###############################		
		argI_i = I_i + 1d0/2d0 * k2_I_i
		argS_i = S_i + 1d0/2d0 * k2_S_i
		argRI_ij = RI_ij + 1d0/2d0 * k2_RI_ij
                argSI_ij = SI_ij + 1d0/2d0 * k2_SI_ij			
		argSS_ij = SS_ij + 1d0/2d0 * k2_SS_ij
		!###############################


		k3_I_i = dt * fI_i(t + (dt/2d0), this, m_sitios, n_arestas, lambda, mu, argI_i, argSI_ij)		
		
		k3_S_i = dt * fS_i(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, argI_i, argS_i, argSI_ij)
		
		k3_RI_ij = dt * fRI_ij(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k3_SI_ij = dt * fSI_ij(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k3_SS_ij = dt * fSS_ij(t + (dt/2d0), this, m_sitios, n_arestas, alp, lambda, mu, argS_i, argSI_ij, argSS_ij)



		!###############################		
		argI_i = I_i + k3_I_i
		argS_i = S_i + k3_S_i
		argRI_ij = RI_ij + k3_RI_ij
                argSI_ij = SI_ij + k3_SI_ij			
		argSS_ij = SS_ij + k3_SS_ij
		!###############################



		k4_I_i = dt * fI_i(t + dt, this, m_sitios, n_arestas, lambda, mu, argI_i, argSI_ij)		
		
		k4_S_i = dt * fS_i(t + dt, this, m_sitios, n_arestas, alp, lambda, argI_i, argS_i, argSI_ij)
		
		k4_RI_ij = dt * fRI_ij(t + dt, this, m_sitios, n_arestas, alp, lambda, mu, argI_i, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k4_SI_ij = dt * fSI_ij(t + dt, this, m_sitios, n_arestas, alp, lambda, mu, argS_i, argRI_ij, argSI_ij, argSS_ij)
		
		k4_SS_ij = dt * fSS_ij(t + dt, this, m_sitios, n_arestas, alp, lambda, mu, argS_i, argSI_ij, argSS_ij)


		I_i = I_i + 1d0/6d0 * (k1_I_i + 2d0 * k2_I_i + 2d0 * k3_I_i + k4_I_i)  		

		S_i = S_i + 1d0/6d0 * (k1_S_i + 2d0 * k2_S_i + 2d0 * k3_S_i + k4_S_i)  		
  				
		RI_ij = RI_ij + 1d0/6d0 * (k1_RI_ij + 2d0 * k2_RI_ij + 2d0 * k3_RI_ij + k4_RI_ij)

		SI_ij = SI_ij + 1d0/6d0 * (k1_SI_ij + 2d0 * k2_SI_ij + 2d0 * k3_SI_ij + k4_SI_ij)		  		

		SS_ij = SS_ij + 1d0/6d0 * (k1_SS_ij + 2d0 * k2_SS_ij + 2d0 * k3_SS_ij + k4_SS_ij)
						
      !#########################SUBROTINA#############################      
          
   end subroutine
   
end module

!#######################
