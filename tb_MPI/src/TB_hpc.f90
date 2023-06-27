! *********************************************************************************
!                             PROGRAMA TB_hpc v0.2
! Curso de HPC II
! Professor: Dr. Rodrigo G. Amorim
! Autor: Mateus Corradini Lopes
!
! O programa tem fins didáticos, para treinar o paralelismo e para intruduzir a
! ideia do método semi-empírico TB.
! Calcula a estrutura de banda de uma cadeia 1D com NAT atomos, todos da
! mesma espécie, com dois orbitais para cada atomo SPx. O sistema esta sujeito a
! condições de contorno periodicas, o último átomo visualiza o primeiro, e 
! aproximação de primeiros vizinhos.
! O programa le o arquivo de entrada input.dat e gera autovalores.dat.
! *********************************************************************************

program TB_hpc
use Sort
implicit none
include "mpif.h"


! -------------------- Declaração de Var. Array --------------------------------

   ! Do método de TB

     integer                  :: i, j, k, nk, ik, flp, dist
     integer                  :: nat, ndim
     real*8                   :: alfa(2), beta(2,2), a, qi,ac,tau
     real*8, allocatable      :: kvec(:), w(:,:)
     complex*8,allocatable    :: H(:,:,:), uk(:,:,:)

  ! Da API MPI

    integer                   :: myid,nprocs,root
	integer                   :: ierr
    integer status(MPI_STATUS_SIZE)


  ! Arquivos

      open(4,file = 'input.dat',status='unknown')
      open(9,file = 'autovalores.dat',status='unknown')
      !open(12,file = 'autoestados.dat', status='unknown')

!------------------ Abertura da sessão paralela  ---------------------------------

    	call MPI_INIT(ierr)
	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
	call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)


!------------------ Leitura dos dados de entrada ---------------------------------
! Feita apenas pelo root.

    if(myid == 0) then

        read(4,*) nat
        read(4,*) alfa(1), alfa(2)
        read(4,*) beta(1,1), beta(2,2), beta(1,2)
        beta(2,1) = beta(1,2)
        read(4,*) a
        read(4,*) nk


        print*
        write(*,'(a,i10)') ' Numero de Atomos na Cadeia(NAT):  ', nat
        write(*,'(a,f10.4,f10.4)') '  On-site Energys (ALFA(1) e ALFA(2)):  ', alfa(1), alfa(2)
        write(*,'(a,f10.4,f10.4,f10.4)') '  Hoppings de primeiros vizinhos(beta):  ', beta(1,1), beta(2,2), beta(1,2)
        write(*,'(a,f10.4)') '  Parâmetro de rede (a): ', a
        write(*,'(a,i10)') '  Número de pontos k (nk): ', nk
        print*

        ndim = 2*nat ! 2 Orbitais/átomo
        flp = nat    ! estado central

    ! numero de matrizes que cada processador vai pegar
    ! + adaptação da grade para divisões não per-
    ! feitas da grade pros processadores.

        dist = int(nk/nprocs)
        nk = nprocs*dist
        write(*,*) " Nova grade nk adaptada: ", nk

    endif

    write(*,*) "Processador: ",myid," esta participando."

!----------------- Alocacao de Mem. + Construção do Hamiltoniano -------------------

    ! Todos os processos precisam saber nk, ndim e dist
    ! e alocar as os arrays necessarios.

    call MPI_BCAST(nk,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ndim,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dist,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

    allocate(H(nk,ndim,ndim), w(nk,ndim),kvec(nk),uk(nk,ndim,ndim))

    if(myid == 0) then

        ! Criando a grade de pontos k

        call gridk(nk,kvec,a)


        ! Chamando a subrotina que constroi o Hamiltoniando de Tight-Binding

        call buildH(nk, H, ndim, nat, alfa, beta, kvec,a)


        ! Distribuindo as matrizes
        do i = 1, nprocs-1
            do j = (i*dist) + 1 , (i+1)*dist
                call MPI_SEND(H(j,1:ndim,1:ndim),ndim*ndim,MPI_COMPLEX,&
                               i,0,MPI_COMM_WORLD,ierr)
            enddo
        enddo

    else

        ! Recebendo as matrizes nos processadores
        do j = (myid*dist)+1, (myid+1)*dist
                call MPI_RECV(H(j,1:ndim,1:ndim),nk*ndim*ndim,MPI_COMPLEX,&
                            0,0,MPI_COMM_WORLD,status,ierr)
        enddo

    endif

    ! Barreira para garantir que todos tenham lidado com suas
    ! matrizes.

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! Rotina da diagonalização e iteração QR para os autovalores
    ! sendo efetuada na sua respectiva faixa em cada processador


    do ik = (myid*dist)+1, (myid+1)*dist

        call qrit(H,w,ndim,nk,ik,uk,myid)

    enddo


    if(myid /= 0) then

            do j = (myid*dist)+1 , (myid+1)*dist
                call MPI_SEND(w(j,1:ndim),ndim,MPI_DOUBLE_PRECISION,&
                               0,1,MPI_COMM_WORLD,ierr)
            enddo

    else

        do i = 1, nprocs-1
            do j = (i*dist)+1, (i+1)*dist
                call MPI_RECV(w(j,1:ndim),nk*ndim,MPI_DOUBLE_PRECISION,&
                                i,1,MPI_COMM_WORLD,status,ierr)
            enddo
        enddo

    endif

    ! Barreira para garantir que todos os autovalores
    ! foram devidamente recebidos no root antes de
    ! finalizazr.

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(myid == 0) then
     !Escreve os autovalores no arquivo

        do i = 1, nk
            call Shell_Sort(w(i,:))
            write(9,'(f10.4,9f10.4)') kvec(i), (w(i,j), j=1 ,2), (w(i,j), j = flp - 2, flp + 2),&
                                        (w(i,j), j = ndim - 1, ndim)
            !write(12,*) kvec(i), (uk(i,j,:), j=1 ,2), (uk(i,j,:), j = flp - 2, flp + 2), (uk(i,j,:), j = ndim - 1, ndim)
            !write(9,'(f10.4,9f10.4)') kvec(i), w(i,:)
        enddo

    endif

    call MPI_FINALIZE(ierr)

    close(4)
    close(9)  
    close(11)

end program TB_hpc

!---------------------------------SUBROTINAS----------------------------------------

    subroutine BuildH(nk, H, ndim, nat, alfa, beta ,kvec,a)
    implicit none
    
    integer         :: i, io1, io2, jo1, jo2, ko1, ko2, ik
    integer, intent(IN) :: ndim, nat , nk
    complex*8,intent(OUT)  :: H(nk,ndim,ndim)
    real*8,intent(IN)   :: alfa(2), beta(2,2) , kvec(nk), a
    complex*8          :: ek

    do ik =1, nk

    ek = cmplx(cos(kvec(ik)*a),sin(kvec(ik)*a))

    H(ik,:,:) = cmplx(0.0d0,0.0d0)

        do i = 1, nat-1  ! Loop sobre os atomos

        ! Na construcao da hamiltoniana, cada atomo i  irá interagir com seus
        ! Vizinhos a direita, de rotulo i + 1. 
        ! Para cada atomo i duas linhas da Hamiltoniana sao construidas, 
        ! uma para cada orbital.

            io1 = (i - 1)*2 + 1 ! Orbital 1, atomo i
            io2 = io1 + 1       ! Orbital 2, atomo i
            jo1 = 2*i + 1       ! Orbital 1, atomo j (primeiro vizinho de i)
            jo2 = jo1 + 1       ! Orbital 2, atomo j (primeiro vizinho de i)

        
            H(ik,io1,io1) = alfa(1) + ek*beta(1,1)   ! Onsite energy principal, orbital 1
            H(ik,io2,io2) = alfa(2) + ek*beta(2,2)   ! Onsite energy principal, orbital 2
            H(ik,io1,jo1) = ek*beta(1,1)  ! Hopping <io1|H|jo1>
            H(ik,io1,jo2) = ek*beta(1,2)  ! Hopping <io1|H|jo2>
            H(ik,io2,jo1) = ek*beta(2,1)  ! Hopping <io2|H|jo1>
            H(ik,io2,jo2) = ek*beta(2,2)  ! Hopping <io2|H|jo2>
         


        ! Impondo hermiticidade aos hoppings:

            H(ik,jo1,io1) = conjg(H(ik,io1,jo1))
            H(ik,jo1,io2) = conjg(H(ik,io2,jo1))
            H(ik,jo2,io1) = conjg(H(ik,io1,jo2))
            H(ik,jo2,io2) = conjg(H(ik,io2,jo2))
        

        enddo

    ! Definindo os elementos de H para i = NAT e impondo condicoes de contorno 
    ! Periodicas: o átomo NAT interage com o átomo i=1 .
 
        io1 = (nat - 1)*2 + 1
        io2 = io1 + 1
        jo1 = 1    
        jo2 = 2


        H(ik,io1,io1) = alfa(1) + ek*beta(1,1)
        H(ik,io2,io2) = alfa(2) + ek*beta(2,2)
        H(ik,io1,jo1) = ek*beta(1,1)
        H(ik,io1,jo2) = ek*beta(1,2)
        H(ik,io2,jo1) = ek*beta(2,1)
        H(ik,io2,jo2) = ek*beta(2,2)

    ! Hermiticidade:

        H(ik,jo1,io1) = conjg(H(ik,io1,jo1))
        H(ik,jo1,io2) = conjg(H(ik,io2,jo1))
        H(ik,jo2,io1) = conjg(H(ik,io1,jo2))
        H(ik,jo2,io2) = conjg(H(ik,io2,jo2))

   !write(*,*) 'H(',ik,')=', H(ik,:,:)

    enddo

    end subroutine buildH

    !------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------



    subroutine qrit(H,w,ndim,nk,ik,uk,myid)
    ! Variáveis externas

    integer,intent(IN)          :: ndim,nk,ik, myid
    real*8,intent(INOUT)        :: w(nk,ndim)
    complex*8,intent(INOUT)     :: H(nk,ndim,ndim), uk(nk,ndim,ndim)

    ! Variáveis internas
    complex*8                  :: u(ndim,ndim), q(ndim,ndim), R(ndim,ndim)
    complex*8                  :: acproj(ndim), proj(ndim)
    integer                    :: i, j, it
    real*8                     :: epsilon, lim
    complex*8                  :: t1(ndim), t2(ndim), a(ndim,ndim)

    epsilon = 1.0d0
    lim = 0.0001
    it = 0

    ! uk guardará os autoestados
!     uk(ik,:,:) = 0.0d0
!     do i=1, ndim
!         uk(ik,i,i) = dcmplx(1.0d0,0.0d0)
!     enddo

!----------------- Diagonalização da Matriz H -------------------

    ! Início das iterações QR

    do while(epsilon > lim)

        ! Reiniciando os valores da iteração QR
        u(:,:) = 0.0d0
        q(:,:) = 0.0d0
        R(:,:) = 0.0d0

        ! -----------Início do Gram-Schmidt modificado-----------
        !  A primeira e segunda colunas sao feitas fora do loop

        u(:,1) = H(ik,:,1)
        q(:,1) = u(:,1)/sqrt(dot_product(u(:,1),u(:,1)))

        proj = (dot_product(u(:,1),H(ik,:,2))/dot_product(u(:,1),u(:,1)))*u(:,1)
        u(:,2) = H(ik,:,2) - proj

        q(:,2) = u(:,2)/sqrt(dot_product(u(:,2),u(:,2)))

        do i = 3, ndim

            acproj = 0.0d0

            proj = (dot_product(u(:,1),H(ik,:,i))/dot_product(u(:,1),u(:,1)))*u(:,1)
            u(:,i) = H(ik,:,i) - proj

            do j = 2, i - 1

                proj = 0.0d0
                call proj_ua(u, u, ndim, proj,i,j)
                u(:,i) = u(:,i) - proj

            enddo

            q(:,i) = u(:,i)/sqrt(dot_product(u(:,i),u(:,i)))

        enddo

        ! Atualiza a matriz dos autoestados

        !call multMatrix(uk(ik,:,:),uk(ik,:,:),q,ndim,ndim,ndim,ndim)

        ! Montando a matriz A

        do i = 1, ndim

            a(:,i) = 0.0d0

            do j = 1, i
                    a(:,i) = a(:,i) + dot_product(q(:,j),H(ik,:,i))*q(:,j)
            enddo

        enddo

        ! Montando a matriz R

        do i = 1, ndim
            do j = i, ndim
                    R(i,j) = dot_product(q(:,i),a(:,j))
            enddo
        enddo

    call multMatrix(H(ik,:,:),R,q,ndim,ndim,ndim,ndim)

    ! Buscando maior valor absoluto abaixo da diagonal principal

        epsilon = 0.0d0

        do i = 2, ndim
            do j = 1, i-1
                if( cabs(H(ik,i,j)) > epsilon) then
                    epsilon = cabs(H(ik,i,j))
                endif
            enddo
        enddo


    it = it + 1
    if(it >= 500000) exit

    end do


    ! Matriz H convergiu, cria o vetor dos autovalores w

    write(*, *) myid, it, epsilon

    do i = 1, ndim
        w(ik,i) = H(ik,i,i)
        !write(*,*) ik, H(ik,i,i)
    enddo

    end subroutine qrit


    !------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------


    subroutine proj_ua(u,H,ndim,proj,i,j)
    implicit none
    
    integer,intent(IN)         ::  ndim, i , j
    complex*8,intent(IN)       ::  u(ndim,ndim)
    complex*8,intent(INOUT)    ::  proj(ndim)
    complex*8,intent(IN)       ::  H(ndim,ndim)
    complex*8                  ::  t1(ndim),t2(ndim)
    complex*8                  ::  sigma
    
    sigma = dot_product(u(:,j),H(:,i))/dot_product(u(:,j),u(:,j))
    proj = sigma * u(:,j)


    end subroutine proj_ua

    !------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------

    subroutine multMatrix(Mres, M1, M2, m1d1, m1d2, m2d1, m2d2)
    implicit none

    integer,intent(IN)  :: m1d1, m1d2, m2d1, m2d2
    complex*8,intent(IN)   :: M1(m1d1,m1d2), M2(m2d1,m2d2)
    complex*8,intent(OUT)  :: Mres(m1d1,m2d2)
    integer     :: i, j, k
    Mres = 0.0d0

    do i=1, m2d2
       do j = 1, m1d2
          do k = 1, m1d1
            Mres(i,j) = Mres(i,j) + (M1(i,k)*M2(k,j))
          end do
       end do
    end do

    end subroutine multMatrix

    !------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------

    subroutine gridk(nk,kvec,a)
    implicit none

    ! Variáveis externas

    integer,intent(IN)      :: nk
    real*8,intent(INOUT)    :: kvec(nk)
    real*8,intent(IN)       :: a

    ! Variáveis internas

    integer                 :: i
    real*8                  :: pi, step, full

    ! Dividindo o range k \in [-pi/a,pi/a]

    pi = 4.d0*atan(1.0d0)
    step = ((2.0d0 * pi)/a)/nk

    do i = 1, nk
        kvec(i) = -(pi/a) + (i-1)*step
    enddo

    end subroutine gridk
