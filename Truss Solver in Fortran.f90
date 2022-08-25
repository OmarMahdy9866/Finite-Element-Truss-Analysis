 PROGRAM TrussSolver
    !MS$DEBUG
    IMPLICIT NONE
    INTERFACE
    SUBROUTINE GeneralStiffness (E,A, GDof, NodesCoordinates, numberElements,numberNodes,elementNodes, Stiffness,ADof) 
    ! Generates the main stiffness matrix of the whole system
    !declaration of function arguments
    real, INTENT(IN) :: E,A
    INTEGER, INTENT(IN) :: numberElements,numberNodes,ADof
    REAL, dimension(1:numberNodes,1:2), INTENT(IN) :: NodesCoordinates
    INTEGER, INTENT(IN) :: GDof
    INTEGER, dimension(1:numberElements,1:2), INTENT(IN) :: elementNodes
    REAL, dimension(1:GDof,1:GDof), INTENT(OUT) :: Stiffness
    END SUBROUTINE GeneralStiffness

    SUBROUTINE  REDUCED_SYSTEM(ADof,GDof_Array,GDof,PDof_ARRAY,PDof,stiffness,stiffness_red,force,force_red,ADof_array) !Generates the reduced system equations
    integer, intent(in) :: GDof,PDof, ADof
    integer, dimension (1:Gdof), intent(in):: GDof_Array
    integer, dimension (1:Pdof), intent(in):: PDof_Array
    real, dimension(1:GDof,1:Gdof), intent(in):: stiffness
    real, dimension(1:GDof), intent(IN) :: force
    real, dimension(1:ADof), intent(OUT) :: force_red
    real, dimension(1:ADof,1:Adof), intent(INOUT):: stiffness_red
    integer,dimension(1:ADof),intent(inout)::ADof_array
    !declaration of function arguments
    END SUBROUTINE REDUCED_SYSTEM
    
    subroutine inverse(a,c,n)
    implicit none 
     real,dimension(n,n), intent(inout):: a
     real, dimension(n,n), intent(OUT) :: c
    integer, intent(in) :: n
    end subroutine inverse
    
    SUBROUTINE Stress(numberNodes, GDof, numberElements,elementNodes,xx,yy,displacement,E, Stresses,Axial_Forces,A)
    implicit none
    integer, intent(in) :: numberelements, GDof, numberNodes
    real, intent(in)::E,A
    real, dimension(1:GDof), intent(in):: displacement
    integer,dimension(1:numberElements,1:2), intent(in):: elementNodes
    real, dimension(1:numberNodes), intent(in) :: xx
    real, dimension(1:numberNodes), intent(in) :: yy
    real, dimension (1:numberelements), intent(inout) :: Stresses
    real, dimension (1:numberelements), intent(inout) :: Axial_Forces
    end subroutine stress
    
    SUBROUTINE SHOW_OUTPUT(GDof,Stresses,Stiffness,displacement, Axial_Forces, numberElements)
    implicit none
    REAL, intent(IN) :: stiffness(GDof,GDof)
    real, intent(in) :: displacement(GDof)
    real, intent(in) :: Axial_Forces(numberElements)
    real, intent(in) :: Stresses(numberElements)
    integer, intent(in) :: numberElements
    integer, intent(in) :: GDof
    END SUBROUTINE SHOW_OUTPUT
    
   END INTERFACE
    !___________________________________________________________________________________________
    !Declaration of Main program variables
    integer :: GDof, n, m, j, t,r,w,i,numberElements,numberNodes, forcesNumber, PDof,ADof
    real, dimension(:), allocatable :: Stresses
    real, dimension(:,:),allocatable::Reactions
    real, dimension(:,:), allocatable ::stiffness   
    real, dimension(:), allocatable :: xx
    real, dimension (:,:), allocatable :: NodesCoordinates
    real, dimension(:), allocatable :: yy
    real, dimension(1:4,1:4) :: K_element=0
      real ,dimension(:,:), allocatable :: Stiffness_Red
    real, dimension(:), allocatable :: Force_Red
    real, dimension(:,:), allocatable :: InvOfStiffRed
    real, dimension(:), allocatable :: Force
    real,dimension(:,:), allocatable :: forceSystem
    real, dimension(:), allocatable :: displacement
     real, dimension(:), allocatable :: displacement_Red
    integer, dimension(:,:),allocatable :: elementNodes !Allocated
    integer, dimension(:), allocatable :: GDof_Array    !Allocated
    integer, dimension(:), allocatable :: PDof_Array    !Allocated
    integer, dimension(:), allocatable :: ADof_Array    !Allocated
    real, dimension(:), allocatable :: Axial_Forces
    real :: E,A,EA,C,S,element_length,deltaX,deltaY,K,t1,t2
    call cpu_time(t1)
     open (unit=9866, file='TRUSS.txt', status='old', action='read')
    read(9866, *) E
    read(9866, *) A
    read(9866, *) numberElements
    allocate(elementNodes(1:numberElements,1:2))
    do i=1,numberElements
        read(9866,*) elementNodes(i,:)
    end do
    numberNodes=maxval(elementNodes)
    allocate(NodesCoordinates(1:numberNodes,1:2))
    do i=1,numberNodes
        read(9866,*) nodesCoordinates(i,:)
    end do
    read(9866,*) forcesNumber ! To be decomposited
    allocate(forceSystem(1:forcesNumber,1:2))
    do i=1,forcesNumber
        read(9866,*) forceSystem(i,:) !This is a combination of the magnitude..
    end do                            !of forces and the DOF associated with it
    read(9866,*) PDof
    allocate(PDof_Array(1:PDof))
    do i=1,PDof
        read(9866,*) PDof_Array(i)
    end do
    
    GDof = 2*maxval(elementNodes)
    ADof=GDof-PDof
    
    allocate(GDof_Array(1:GDof))
    allocate(stiffness_red(1:ADof,1:ADof))
    allocate(displacement_Red(1:ADof))
    allocate(Force_Red(1:ADof))
    allocate(stiffness(1:GDof,1:GDof))
    allocate(force(1:Gdof))
    allocate(displacement(1:Gdof))
    allocate(InvOfStiffRed(1:ADof,1:ADof))
    allocate(Stresses(numberElements))
    allocate(Axial_Forces(numberElements))
    allocate(adof_array(adof))
    allocate(reactions(PDof,2))
      xx=NodesCoordinates(:,1)
      yy=NodesCoordinates(:,2)
    
    Stresses=0
    stiffness=0
    force=0
    displacement=0
    stiffness_red=0.00
    force_Red=0
    
    do i=1,forcesNumber
    n=int(forceSystem(i,1))
        force(n)=forceSystem(i,2)
    end do
    do i=1,GDof
        GDof_Array(i)=i
    end do

    CALL GeneralStiffness (E,A, GDof, NodesCoordinates, numberElements,numberNodes,elementNodes, Stiffness,ADof)
    
    CALL REDUCED_SYSTEM(ADof, GDof_Array,GDof,PDof_ARRAY,PDof,stiffness,stiffness_red,force,force_red,ADof_ARRAY)
   
    CALL inverse(stiffness_red,invofstiffred,adof)
    
    displacement_red=matmul(invofstiffred,force_red)
    
    do i=1,GDof
        do j=1,ADof
            m=ADof_ARRAY(j)
            if (i==m) then
            displacement(i) = displacement_red(j)
           end if
           end do
           end do
    
    CALL Stress(numberNodes, GDof, numberElements,elementNodes,xx,yy,displacement,E, Stresses,Axial_Forces,A)
    
    force=matmul(stiffness,displacement)
    
     do i=1,GDof
        do j=1,PDof
        m=PDof_ARRAY(j)
        if (i==m) then
            reactions(j,1)=m
            reactions(j,2)=force(i)
        end if
         end do
        end do
    do i=1,numberElements
        write(*,*)Axial_Forces(i)
    end do
   ! CALL ExtractResults(Stresses,Stiffness,Displacement,Force,Displacement_Red,Stiffness_Red,Force_Red,invofstiffred,GDof,ADof)
    CALL SHOW_OUTPUT(GDof,Stresses,Stiffness,displacement, Axial_Forces, numberElements)
    call cpu_time(t2)
    write(*,*)"Program starts ", t1
    write(*,*)"Program ends ", t2
    write(*,*)"Excution time =",t2-t1
    
 end program TRUSSSOLVER
    !_______________________________________________________________________________________________________________PROGRAM ENDS HERE
    
    SUBROUTINE GeneralStiffness (E,A, GDof, NodesCoordinates, numberElements, numberNodes, elementNodes, Stiffness,ADof)
    real, INTENT(IN) :: E,A
    INTEGER, INTENT(IN) :: GDof,ADof
    REAL, dimension(1:numberNodes,1:2), INTENT(IN) :: NodesCoordinates
    INTEGER, INTENT(IN) :: numberElements,numberNodes
    INTEGER, dimension(1:numberElements,1:2), INTENT(IN) :: elementNodes
    REAL, dimension(1:GDof,1:GDof), INTENT(OUT) :: Stiffness
    
    !__________________________________________________________________
    
    integer :: z, n, m, j, t,r,w,i
    real, dimension(1:NumberNodes) :: xx
    real, dimension(1:NumberNodes) :: yy 
    real, dimension(1:4,1:4)::K_element
    integer, dimension(1,1:4) :: elementDof
    integer, dimension(1,1:2) :: indexx
    real :: C,S,element_length,deltaX,deltaY,K
    stiffness=0

    xx=NodesCoordinates(:,1)
    yy=NodesCoordinates(:,2)

!Input Check Step
open(3, file = 'INPUTSCHECK.txt', status = 'replace')
write(3,*) "NodesCoordinates = "
do i=1,numberNodes  
    write(3,*)(NODESCOORDINATES(i,j),j=1,2)
end do
write(3,*)"E = ",E
write(3,*)"A = ",A
write(3,*)"GDof = ",GDOF
write(3,*)"ADof=",ADof
write(3,*)"Number of elements = ",numberElements
write(3,*)"numberNodes = ",numberNodes
write(3,*)"Elementnodes = "

do i=1,numberelements
write(3,*)(elementNodes(i,j),j=1,2)
end do
do i=1,numberElements
        do j=1,2    ! This is, only, to fill in the indexx array
         indexx(1,j)=elementNodes(i,j) !Connectivity Array of each element
        end do
    elementDof=reshape((/indexx(1,1)*2-1, indexx(1,1)*2, indexx(1,2)*2-1, indexx(1,2)*2/),(/1,4/))
    deltaX=xx(indexx(1,2))-xx(indexx(1,1))
    deltaY=yy(indexx(1,2))-yy(indexx(1,1))
    element_length=(deltaX**2+deltaY**2)**0.5
    C=deltaX/element_length
    S=deltaY/element_length
    K=E*A/element_length
    K_ELEMENT=K*reshape((/C**2, C*S, -C*C, -C*S, C*S, S*S, -C*S, -S*S, -C*C, -C*S, C*C, C*S, -C*S, -S*S, C*S, S*S/),(/4,4/)) !Symmetrix matix
    
    !!! Assemblage step
        do z=1,4
             n=elementDof(1,z)
             do w=1,4
                m=elementDof(1,w)
                stiffness(n,m)=K_ELEMENT(z,w)+stiffness(n,m)
            end do
         end do
    end do
    !__________________________________________________________________________________________________________ !!!!! Assemblage Done
open (unit=1, file='Stiffness.xls', status='replace')
do j=1,GDof
write(1,*)(stiffness(j,r),r=1,Gdof)
end do

END SUBROUTINE GeneralStiffness
    
    SUBROUTINE REDUCED_SYSTEM(ADof,GDof_Array,GDof,PDof_ARRAY,PDof,stiffness,stiffness_red,force,force_red,ADof_ARRAY)
    integer, intent(in) :: GDof,PDof, ADof
    integer, dimension (1:Gdof), intent(in):: GDof_Array
    integer, dimension (1:Pdof), intent(in):: PDof_Array
    real, dimension(1:GDof,1:Gdof), intent(in):: stiffness
    real, dimension(1:GDof), intent(IN) :: force
    real, dimension(1:ADof), intent(OUT) :: force_red
    real, dimension(1:ADof,1:Adof), intent(INOUT):: stiffness_red
    integer,dimension (1:Adof), intent(INOUT) :: ADof_Array
    !___________________________________________________________
    
    integer :: i,j,n,m
    
    
    ADof_Array=pack(GDof_ARRAY, [(all(PDof_ARRAY /= GDof_ARRAY(i)), i=1, size(GDof_ARRAY))])
    stiffness_red=0
    do i=1,ADof
        n=ADOF_ARRAY(i)
            do j=1,ADof
                m=ADof_ARRAY(j)
                Stiffness_Red(i,j)=Stiffness_Red(i,j)+stiffness(n,m)
            end do
    end do
    
    do i=1,ADof
         n=ADof_ARRAY(i)
         force_red(i)=force(n)
    end do

    end Subroutine REDUCED_SYSTEM
    
    
    subroutine Inverse(a,c,n)
    implicit none 
    real,dimension(n,n), intent(inout):: a
     real, dimension(n,n), intent(OUT) :: c
    integer, intent(in) :: n
     real :: L(n,n), U(n,n), b(n), d(n), x(n)
    real :: coeff
    integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
    end subroutine inverse
    
    SUBROUTINE Stress(numberNodes, GDof, numberElements,elementNodes,xx,yy,displacement,E,Stresses,Axial_Forces,A)
    implicit none
    integer, intent(in) :: numberelements, GDof, numberNodes
    real, intent(in)::E,A
    real, dimension(1:GDof), intent(in):: displacement
    integer,dimension(1:numberElements,1:2), intent(in):: elementNodes
    real, dimension(1:numberNodes), intent(in) :: xx
    real, dimension(1:numberNodes), intent(in) :: yy
    real, dimension (numberelements), intent(inout) :: Stresses
    real, dimension (1:numberelements), intent(inout) :: Axial_Forces
    !_________________________________________________________________
    integer::i,j,k
    real:: deltaX, deltaY, element_length, C, S
    real,dimension(4)::trig
    integer, dimension(1,1:4) :: elementDof
    integer, dimension(1,1:2) :: indexx
    real,dimension(4) :: disp
    
    do i=1,numberElements
        do j=1,2    ! This is, only, to fill in the indexx array
         indexx(1,j)=elementNodes(i,j) !Connectivity Array of each element
        end do
    elementDof=reshape((/indexx(1,1)*2-1, indexx(1,1)*2, indexx(1,2)*2-1, indexx(1,2)*2/),(/1,4/))
    deltaX=xx(indexx(1,2))-xx(indexx(1,1))
    deltaY=yy(indexx(1,2))-yy(indexx(1,1))
    element_length=(deltaX**2+deltaY**2)**0.5
    C=deltaX/element_length
    S=deltaY/element_length
    Trig=(/-C,-S,C,S/)
    disp=displacement(elementDof(1,1:4))
    !stresses(i)=E/element_length*(-C*elementDof(1,1)-S*elementDof(1,2)+C*elementDof(1,3)+S*elementDof(1,4))
    stresses(i)=E/element_length*dot_product(disp,trig)
    Axial_Forces(i)=Stresses(i)*A
    end do
    
    end subroutine
    

SUBROUTINE SHOW_OUTPUT(GDof,Stresses,Stiffness,displacement, Axial_Forces, numberElements)
    implicit none
    
      REAL, intent(IN) :: stiffness(GDof,GDof)
    real, intent(in) :: displacement(GDof)
    real, intent(in) :: Axial_Forces(numberElements)
    real, intent(in) :: Stresses(numberElements)
    integer, intent(in) :: numberElements
    integer, intent(in) :: GDof
    
    integer:: i,j
     open (unit=7, file='OUTPUT.xls', status='replace')
     
     write(7,*)"Golbal Stiffness Matrix (K) "
     write(7,*)" "
     do i=1,GDof
     write(7,*)(Stiffness(i,j),j=1,GDof)
     end do
     write(7,*)" "
     write(7,*)"Axial Forces in Each member "
     write(7,*)" "
     do i=1,numberElements
     write(7,*)Axial_Forces(i)
     end do
     write(7,*)" "
     write(7,*)" The Displacement at Each Node "
     write(7,*)" "
     do i=1,GDof
     write(7,*)Displacement(i)
     end do
     write(7,*)" "
     write(7,*)"The Stresses in Each Element "
     write(7,*)" "
     do i=1,numberElements
     write(7,*)Stresses(i)
     end do
end subroutine SHOW_OUTPUT