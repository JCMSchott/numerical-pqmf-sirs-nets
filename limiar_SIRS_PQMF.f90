program main
   use geraRede
   use mod_rndgen
   use mod_tools
   implicit none
   !====================================================================
   type(grafo_PL_UCM) :: rede
   !====================================================================
   integer, parameter :: dp = kind(0.0d0)
   real(dp) :: lamb
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
   !====================================================================
   integer :: seed(10)
   integer :: seed1
   type(rndgen) :: ger_inic
   !====================================================================
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   real(dp) :: grau_max      
   !====================================================================
   character(len=500) :: lamb_vs_Im
   character(len=500) :: lamb_vs_Xi 
   character(len=1000) :: caminho
   character(len=500) :: arquivo_rede_real
   character(len=5) :: tam_char
   character(len=5) :: gama_char
   character(len=5) :: indice
   !====================================================================   
   integer :: ind_am
   logical :: T_vs
   !====================================================================
   integer :: sumDeg2
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   !====================================================================
   integer :: n_it, per_conv
   !====================================================================
   real(dp) :: p0
   character(len=300) :: cwd, resultados
   character(len=1000) :: local
   character(len=100) :: buffer
   !====================================================================   
   real(dp) :: qm, q2m
   !====================================================================   
   integer :: nargus, ind_lamb   
   character(len=3) :: char_ind_lamb
   !====================================================================
   integer, parameter  :: N_iteracoes = 1000
   !====================================================================   
   real(dp) :: a_par
   !====================================================================
   real(dp) :: autovalor_minus, autovalor_plus, autovalor0, autovalor
   !====================================================================
   integer :: i1, i2
   !====================================================================  
   real(dp), allocatable :: Adj(:)
   real(dp), allocatable :: x(:), y(:)
   !====================================================================
   real(dp) :: erro
   real(dp) :: tolet
   real(dp) :: tole
   !====================================================================
   real(dp) :: lambda0, lambda_backup, dlambda0, dlambda
   real(dp) :: lambda_plus, lambda_minus
   real(dp) :: lambda_mean
   real(dp) ::lambda_mean0
   !====================================================================   
   logical :: JaFoiPlus, JaFoiMinus   
   !====================================================================
   integer :: semente2
   type(rndgen) :: geni
   !====================================================================

   !====================================================================
   semente2 = 987652814 
   !====================================================================
   call geni%init(semente2)
   !====================================================================
   
   !====================================================================
   call getcwd(cwd)

   cwd = trim(adjustl(cwd))
   
   resultados = 'Rst_limiar_pQMF'
   resultados = trim(adjustl(resultados))

   call system('mkdir -p '//resultados)
    
   local = trim(adjustl(trim(adjustl(cwd))//"/"//trim(adjustl(resultados))//"/"))

!=======================================================================
   tole = 1.0d-7
   tolet = 1.0d-5
   write(*,*) "================================================================="
   write(*,*) "Tolerancia de ", tole
   write(*,*) "================================================================="
!=======================================================================   
   seed1=947361823
   
   nargus = iargc()

   if(nargus == 5)then
      !#############################
      ! Amostra
      !#############################
      call getarg(1, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) ind_am

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
      ! Alfa
      !#############################
      call getarg(5, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) alp
      write(*,*) "O valor de alp = ", alp
      write(*,*) "================================================================="
   else
      stop "Erro vindo da linha 122 do arquivo 'limiar_SIRS_PQMF.f90'. Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
      write(*,*) "============================================================================================================="
   endif
!=======================================================================
   a_par = (alp * mu)/(alp + 2.0_dp * mu)
!=======================================================================   
   !grau_min = 3
   !tam_rede = 1*10**3
!=======================================================================
   !gama_exp = 2.3_dp
!=======================================================================
!  Indice da semente
   !ind_am = 1   
!=======================================================================
   if(gama_exp < 3.0_dp)then
      grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
   else
      grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
      !grau_max = 2.0_dp * (1.0_dp * tam_rede)**(1d0/(gama_exp-1.d0))
      !call acha_cutoff_rigido(grau_min, gama_exp, tam_rede)
   endif
!=======================================================================
   p0 = 10.0_dp/(1.0_dp * tam_rede) !1.0d-2
!=======================================================================
  write(indice,'(I0)') ind_am   
!=======================================================================
   if(tam_rede == 1000)then
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
   else
      stop 'Escolha um tamanho de rede dentro do catalogo'
      write(*,*) "================================================================="
   endif
!=======================================================================
   write(gama_char,'(f5.2)') gama_exp
!=======================================================================
   call ger_inic%init(seed1)
   i2 = 1
   do i1 = 1, 1000
      if(mod(i1,100) > 0) cycle
      seed(i2)  = ger_inic%int(100000000,999999999)
      write(*,*) i1, seed(i2)
      i2 = i2+1      
   enddo
   
   local = trim(adjustl(trim(adjustl(local))//'/tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
   
   local = trim(adjustl(trim(adjustl(local))//'gam_'//trim(adjustl(gama_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )   

   local = trim(adjustl(trim(adjustl(local))//'ams_'//trim(adjustl(indice))//'/'))
   
   call system('mkdir -p '//trim(adjustl(local)) )


!=======================================================================
!				Inicia grafo
!=======================================================================
 !  arquivo_rede_real='s00088.s838.net.edg'
 !  call rede%RedeReal(arquivo_rede_real, 111)
 !  close(111)
!=======================================================================
  call rede%iniciaGrafo(tam_rede)
  write(*,*) "=========================================================="
!=======================================================================
  call rede%inicia(grau_min, grau_max, gama_exp, seed(ind_am))
  write(*,*) "=========================================================="
!=======================================================================
  indice = trim(adjustl(indice))
!=======================================================================
  indice = trim(adjustl('_'//trim(adjustl(indice)))) 
!=======================================================================
  call rede%liga(seed(ind_am), .False.)  
  write(*,*) "=========================================================="  
!=======================================================================
  call sub_classifica_clusters(rede,.False., 000, 'sem_arquivo.dat')
  write(*,*) "=========================================================="  
!=======================================================================
  write(*,*) "=========================================================="
   write(*,*) ""
   write(*,*) "Tamanho da rede ", rede%nodes, "."
   write(*,*) ""
   write(*,*) "Fracao correspondente aa componente gigante ", 100.0 * comp_gigante/rede%nodes,"%", "."
   write(*,*) ""
   write(*,*) "Expoente da distribuicao de graus da rede ", gama_exp, "."
   write(*,*) ""
   write(*,*) "Grau minimo ", rede%degMin, ".", " Grau maximo ", rede%degMax, "."
   write(*,*) ""
  write(*,*) "==========================================================" 

   arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'.dat'
   call calcula_k_nn(rede, .True., 800, trim(adjustl(arq_1)))
   close(800)

   arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering'))//'.dat'
   call clustering(rede,.True., 800, trim(adjustl(arq_1)))
   close(800)


!#######################################################################
   if(allocated(P_grau)) deallocate(P_grau)
   allocate(P_grau(rede%degMin:rede%degMax))
   P_grau = 0.0_dp

   do i1 = 1, rede%nodes
      P_grau(rede%deg(i1)) = P_grau(rede%deg(i1)) + 1.0_dp
   enddo

   P_grau = P_grau/(1.0_dp * rede%nodes)

   !write(*,*) P_grau

   arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_P_grau_'))//'.dat'

   open(800, file=trim(adjustl(arq_1)), status='unknown')

   do i1 = rede%degMin, rede%degMax
      if(P_grau(i1) == 0.0_dp) cycle
      write(800,*) i1, P_grau(i1)
   enddo
   close(800)
   deallocate(P_grau)
!#######################################################################     
    
   !====================================================================
   write(alp_char2, '(f9.3)') alp
   
   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//'/'))
      
   call system('mkdir -p '//trim(adjustl(local)) )
   
   open(334, file=trim(adjustl(local))//trim(adjustl('lbC_PQMF.dat')), status = 'unknown')
   !====================================================================
   i1 = 0
   !====================================================================
   ! lambda_minus = 0.0d0 corresponde a autolavor = -1
   ! Faz sentido, porque lambda = 0 corresponde ao ponto fixo
   ! I_i = R_i = 0.0 e S_i = 1.0
!======================================================================= 
   call aloca_listas_e_matrizes(rede)
   
   lambda0 = limiarQMF(rede)
   
   dlambda0 = lambda0/10.0d0
   dlambda = dlambda0
     
   lambda_backup = lambda0
   
   JaFoiPlus = .False.
   !====================================================================   
   i1 = 0
   !====================================================================
   lambda_minus = lambda0
   !====================================================================
   ! Procura um lambda_minus apropriado
   do
      i1 = i1 + 1
      lambda_minus = lambda0
!======================================================================= 
      call metodo_potencia(rede, alp, mu, lambda_minus)      
!======================================================================= 
      if( autovalor > 0.0d0 )then
         lambda0 = lambda0 - dlambda
         !write(*,*) "Teve que tirar um pouco do lambda_minus"
      elseif( autovalor == 0.0d0)then
         write(*,*) "Achamos o autovalor nulo por acidente =D ."
         write(*,*) "Lambda_C = ", lambda_minus
         write(334,*) rede%nodes, lambda_minus
         stop
      else
         exit
      endif
      
      if( lambda0 <= 0.0d0 )then
         dlambda = dlambda/2.0d0
         lambda0 = lambda_backup - dlambda
         !write(*,*) "Lambda negativo. Teve que atualizar. dlambda atualizado para ", dlambda
      elseif( lambda0 >= 1.0d0)then
         lambda0 = lambda0 - dlambda0
         dlambda = dlambda0 !/2.0d0
         !write(*,*) "Lambda acima de 1. dlambda atualizado para ", dlambda
      endif
   enddo
   write(*,*) "===================================================================="
   write(*,*) "Achou lambda minus com ", i1, " iteracoes."
!=======================================================================
   write(*,*) "===================================================================="
   write(*,*) "O autovalor associado a lambda_minus = ", lambda_minus, " eh ", autovalor
   write(*,*) "===================================================================="
!======================================================================= 
   autovalor_minus = autovalor
!=======================================================================   
   lambda_plus = 0.10_dp 
!=======================================================================
   i1 = 0
!=======================================================================
   write(*,*) "===================================================================="
   write(*,*) "Vai procurar autovalor_plus"
   write(*,*) "===================================================================="
!=======================================================================  

   lambda0 = limiarQMF(rede)
   
   dlambda0 = lambda0/10.0d0
   dlambda = dlambda0
     
   lambda_backup = lambda0
   
   !====================================================================   
   i1 = 0
   !====================================================================
   lambda_plus = lambda0
   !====================================================================
   ! Procura um lambda_plus apropriado
   do
      i1 = i1 + 1
!======================================================================= 
      call metodo_potencia(rede, alp, mu, lambda_plus)      
!======================================================================= 
      if( autovalor < 0.0d0 )then
         lambda0 = lambda0 + dlambda
         !dlambda = 2.0d0 * dlambda
         !write(*,*) "Dlambda atualizado para ", dlambda
         write(*,*) "Teve que acrescentar um pouco do lambda_plus. Lambda_plus = ", lambda0
      elseif( autovalor == 0.0d0)then
         write(*,*) "Achamos o autovalor nulo por acidente =D ."
         write(*,*) "Lambda_C = ", lambda_plus
         write(334,*) rede%nodes, lambda_plus
         stop
      else
         exit
      endif
      
      !if(lambda0 >= 1.0d0)then
         !lambda0 = lambda0 - dlambda
      !endif
      
      lambda_plus = lambda0
            
   enddo

!=======================================================================   
   autovalor_plus = autovalor
!=======================================================================
   write(*,*) "===================================================================="
   write(*,*) "Achou lambda_plus= ", lambda_plus, " depois de ", i1, " iteracoes."
   
   write(*,*) "===================================================================="
   write(*,*) "Inicializando o metodo da bissecao..."
   write(*,*) "===================================================================="
!=======================================================================
   lambda_mean = (lambda_minus + lambda_plus)/2.0_dp
!=======================================================================
   write(*,*) "O valor de lambda_mean = ", lambda_mean
   write(*,*) "===================================================================="
   i1 = 0   
   do while(i1 < N_iteracoes)
      i1 = i1 + 1
      
!=======================================================================
      call metodo_potencia(rede, alp, mu, lambda_mean)
!=======================================================================
      
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
      if( abs( lambda_plus - lambda_minus )/lambda_mean < tole)then
         write(*,*) "lambda_C = ", lambda_mean, " com erro relativo de ", (100.0d0 * abs(lambda_plus-lambda_minus)/lambda_mean), "%"
         exit
      endif
!=======================================================================      
   enddo
!=======================================================================
   if( i1 < N_iteracoes)then
      write(*,*) "Lambda_C encontrado em ", i1, " iteracoes."
      write(334,*) rede%nodes, lambda_mean
   else
      write(*,*) "Numero de iteracoes excedido."
   endif
   
  contains

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
   
   subroutine metodo_potencia(this, alp1, mu1, lambda1)
      !=================================================================
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp1, mu1, lambda1
      !=================================================================
      integer :: Nmax
      integer :: j1, j2, j3, j4
      real(dp) :: xp, yp
      integer :: ipx, ipy
      type(rndgen) :: geni
      real(dp) :: soma
      integer :: vizim
      real(dp) :: a_l, a_2l
      
      !=================================================================   
      a_l = func_a_l(alp1, lambda1, mu1)
      a_2l = func_a_2l(alp1, lambda1, mu1)      
      !=================================================================
      
      !=================================================================      
      call geni%init(semente2)     
      !=================================================================
      Nmax = 1000
      !=================================================================
	  ! Condicao inicial				
	  !x = 1.0d0
	  !xp = 1.0d0
      !ipx = 1
      
      xp = 0.0d0
      ipx = 1
      do j1 = 1, this%nodes
         x(j1) = geni%rnd() * 1.0d-3/(1.0d0 * this%nodes)
         if( abs(x(j1)) > abs(xp) )then
            xp = x(j1)
            ipx = j1
         endif
      enddo
      x = x/xp
      xp = 1.0d0
      !=================================================================
      ! Apos normalizar o vetor x, a norma ||xp|| = ||x||_{oo}
      ! se torna = 1.0d0.
      !=================================================================
      j1 = 0
      !=================================================================  
!=======================================================================      
iter: do while( j1 < Nmax )	
	     !Aqui comeca o algoritmo
         !==============================================================
         yp = 0.0d0
         ipy = 1
         !==============================================================
         do j2 = 1, this%nodes
            !===========================================================
            soma = 0.0_dp
            do j3 = this%aux(j2), this%aux(j2) + this%deg(j2) - 1
               vizim = this%listAdj(j3)
               soma = soma + x(vizim)
            enddo
            y(j2) = a_l * soma
            !===========================================================
            y(j2) = y(j2) + a_2l * ( 1.0d0 * ( this%degMax - this%deg(j2)) ) * x(j2)
            !===========================================================
            if( abs(y(j2)) > abs(yp) )then
               yp = y(j2)
               !========================================================
               ! Esse ipy sera o proximo ipx.
               !========================================================
               ipy = j2
            endif
            !===========================================================
         enddo
         !write(*,*) "y"
         !write(*,*) y
         !stop
         !==============================================================
         ! A componente y(ipx) representa o ganho multiplicativo
         ! obtido com o produto matricial.
         !==============================================================
         autovalor0 = y(ipx)
		 !==============================================================
		 ! Se yp = 0, para tudo!
	     if(yp == 0.0d0)then
	        xp = 0.0d0
	        ipx = 1
            do j2 = 1, this%nodes   
			   x(j2) = 2.0d0 * geni%rnd()
			   if( abs(x(j2)) > abs(xp) )then
			      xp = x(j2)
			      ipx = j2
			   endif
			enddo
			x = x/xp
			!===========================================================
			! Como o vetor x foi normalizado via || ||_{oo},
			! xp = 1.0d0
			!===========================================================
			xp = 1.0d0
			j1 = 0
			write(*,*) 'Autovalor nulo encontrado. Retomando outro vetor x...'
			cycle iter
         endif
         !==============================================================
		 erro = abs(maxval(x - y/yp))
		 !write(*,*) "O valor de yp eh ", yp
		 
		 !==============================================================
		 x = y/yp
         !==============================================================   
         if(erro < tole)then
            !write(*,*) "Achou autovalor princial = ", autovalor0, " apos ", j1, " iteracoes."
            autovalor = autovalor0 - ( a_2l * (1.0d0 * this%degMax) + mu1 )
            !write(*,*) "autovalor do jacobiano eh ", autovalor               
            exit iter
         endif
         !==============================================================
		 ! O ipy achado apos o produto matricial eh o indice ipx da
		 ! primeira maior componente do vetor x.
		 ! xp neste caso agora = 1.0d0, pois y ja foi normalizado
		 ! na linha anterior ao comando 'if(erro < tole)then'
		 ipx = ipy
		 xp = 1.0d0
         !==============================================================
		 !do j2 = 1, this%nodes
		    !if(abs(x(j2)) > abs(xp))then
		       !xp = x(j2)
		       !ipx = j2
		    !endif
		 !enddo
		 !==============================================================
         j1 = j1 + 1
         !==============================================================
	enddo iter      
      
      
   end subroutine

   function func_a_l(alfa, lambda, mu)
      real(dp), intent(in) :: alfa
      real(dp), intent(in) :: lambda
      real(dp), intent(in) :: mu
      real(dp) :: func_a_l
       
      func_a_l = lambda * ( ( a_par + lambda + mu )/( a_par + 2.0_dp * lambda + mu ) )      
          
   end function
   
   function func_a_2l(alfa, lambda, mu)
      real(dp), intent(in) :: alfa
      real(dp), intent(in) :: lambda
      real(dp), intent(in) :: mu
      real(dp) :: func_a_2l
             
      func_a_2l = (lambda ** 2.0_dp)/( a_par + 2.0_dp * lambda + mu )      
      
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
   
end program
