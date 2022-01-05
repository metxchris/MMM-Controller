!------------------------------------------------------------------------------
!                               Program Description
!------------------------------------------------------------------------------
! mmm_wrapper was designed to be used as a streamlined wrapper around modmmm, as
! an alternative option to testmmm. mmm_wrapper only outputs a single CSV of
! output variables, which is easier to parse than the output produced by
! testmmm.  mmm_wrapper also compiles much quicker than testmmm, which is
! useful for developers who are working on the mmm source files.  While
! mmm_wrapper uses the same input file formatting that testmmm uses, it only
! supports inputs of the first kind; users who wish to use other kinds of
! inputs should still use testmmm.
!
! Additionally, mmm_wrapper serves as a pseudo-replacement for interfacing
! functions to other programming languages, such as the MEX function in Matlab.
! As such, mmm_wrapper assumes the input file is written correctly and offers
! no support for input variable calculations, while also doing minimal error
! checks on values read from the input file.  Users who encounter issues
! running mmm_wrapper should run their input file in testmmm to receive better
! diagnostic feedback.


PROGRAM mmm_wrapper
USE modmmm

IMPLICIT NONE

!------------------------------------------------------------------------------
!                               Local Variables
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(12, 100)

INTEGER, PARAMETER :: &
    hfIn = 34, &  ! Input file handle 
    hfOut = 35    ! Output file handle

INTEGER(8) :: &
    tic, toc, count_rate, count_max  ! Timing variables

! Constants for missing input detection
REAL(R8), PARAMETER :: BADREAL = -1E10_R8
INTEGER, PARAMETER  :: BADINT = -1000000

! Loop iterators
INTEGER :: i, j

!------------------------------------------------------------------------------
!                               Input Controls
!------------------------------------------------------------------------------
INTEGER :: &
    input_kind = BADINT, &  ! Only 1st kind is supported
    npoints = BADINT        ! Number of radial points

!------------------------------------------------------------------------------
!                               Input Variables
!------------------------------------------------------------------------------
REAL(R8), DIMENSION(MMM_NMODE) :: &
    cmodel = BADREAL  ! Internal model weights

REAL(R8) :: &
    cswitch(MAXNOPT, MMM_NMODE)  ! Holds real options for each model

INTEGER :: &
    lswitch(MAXNOPT, MMM_NMODE)  ! Holds integer options for each model

! Real options for each model, to be passed to cswitch (Table 2 in documentation)
REAL(R8), DIMENSION(MAXNOPT) :: & 
    cW20 = BADREAL, &
    cDBM = BADREAL, &
    cETG = BADREAL, &
    cMTM = BADREAL, &
    cETGM = BADREAL 

! Integer options for each model, to be passed to lswitch (Table 3 in documentation)
INTEGER, DIMENSION(MAXNOPT) :: &
    lW20 = BADINT, &
    lDBM = BADINT, &
    lETG = BADINT, &
    lMTM = BADINT, &
    lETGM = BADINT

INTEGER :: &
    lprint = 0  ! Verbose level

REAL(R8), ALLOCATABLE, DIMENSION(:) :: &
    rmin,   &  ! Half-width of the flux surface [m]
    rmaj,   &  ! Major radius to geometric center of the flux surface [m]
    elong,  &  ! Local elongation of flux surface
    ne,     &  ! Electron density [m^-3]
    nh,     &  ! Thermal hydrogenic ion density [m^-3]
    nz,     &  ! Impurity ion density [m^-3]
    nf,     &  ! Density from fast (non-thermal) ions [m^-3]
    zeff,   &  ! Mean Charge
    te,     &  ! Electron temperature [keV]
    ti,     &  ! Temperature of thermal ions [keV]
    q,      &  ! Magnetic q-value, Safety Factor
    btor,   &  ! R0 Bz / R [Tesla]
    zimp,   &  ! Mean charge of impurities
    aimp,   &  ! Mean atomic mass of impurities
    ahyd,   &  ! Mean atomic mass of hydrogen ions
    aimass, &  ! Mean atomic mass of thermal ions
    wexbs,  &  ! ExB shearing rate [s^-1]
    vtor,   &  ! Toroidal velocity [m/s]
    vpol,   &  ! Poloidal velocity [m/s]
    vpar,   &  ! Parallel velocity [m/s]
    gne,    &  ! -R (dne / dr) / ne
    gni,    &  ! -R (dni / dr) / ni
    gnh,    &  ! -R (dnh / dr) / nh
    gnz,    &  ! -R (dnz / dr) / nz
    gte,    &  ! -R (dTe / dr) / Te
    gti,    &  ! -R (dTi / dr) / Ti
    gq,     &  !  R (dq / dr) / q
    gvtor,  &  !  R (dvtor / dr) / vtor
    gvpol,  &  !  R (dvpol / dr) / vpol
    gvpar      !  R (dvpar / dr) / vpar

!------------------------------------------------------------------------------
!                               Output Variables
!------------------------------------------------------------------------------
REAL(R8), ALLOCATABLE, DIMENSION(:) :: &
    xti,        &  ! Effective ion thermal diffusivity [m^2/s]
    xdi,        &  ! Effective hydrogenic ion diffusivity [m^2/s]
    xte,        &  ! Effective electron thermal diffusivity [m^2/s]
    xdz,        &  ! Impurity ion diffusivity from the Weiland model [m^2/s]
    xvt,        &  ! Toroidal momentum transport from the Weiland model [m^2/s]
    xvp,        &  ! Poloidal momentum transport from the Weiland model [m^2/s]
    xtiW20,     &  ! Ion thermal diffusivity from the Weiland model [m^2/s]
    xdiW20,     &  ! Particle diffusivity from the Weiland model [m^2/s]
    xteW20,     &  ! Electron thermal diffusivity from the Weiland model [m^2/s]
    xtiDBM,     &  ! Ion thermal diffusivity from the DRIBM model [m^2/s]
    xdiDBM,     &  ! Hydrogenic ion diffusivity from the DRIBM model [m^2/s]
    xteDBM,     &  ! Electron thermal diffusivity from the DRIBM model [m^2/s]
    xteETG,     &  ! Electron thermal diffusivity from the Horton ETG model [m^2/s]
    xteMTM,     &  ! Electron thermal diffusivity from MTM model [m^2/s]
    xteETGM,    &  ! Electron thermal diffusivity from the ETGM model [m^2/s]
    xdiETGM,    &  ! Hydrogenic ion diffusivity from the ETGM model [m^2/s]
    gammaDBM,   &  ! Growth rate of the most unstable DRIBM mode [s^-1]
    omegaDBM,   &  ! Frequency of the most unstable DRIBM mode [s^-1]
    gammaMTM,   &  ! Growth rate of the most unstable MTM mode [s^-1]
    omegaMTM,   &  ! Frequency of the most unstable MTM mode [s^-1]
    gammaETGM,  &  ! Growth rate of the most unstable ETGM mode [s^-1]
    omegaETGM,  &  ! Frequency of the most unstable ETGM mode [s^-1]
    dbsqprf        ! |(delta B)/B|^2 profile

REAL(R8), ALLOCATABLE, DIMENSION(:, :) :: &
    gammaW20, &  ! Growth rate of the most unstable ion, electron mode in Weiland for +, - directions
    omegaW20     ! Frequency of the most unstable ion, electron mode in Weiland for +, - directions

REAL(R8), ALLOCATABLE, DIMENSION(:, :) :: &
    vflux, &  ! Total flux for each particle type [W/m^2]
    vconv     ! Convective velocities [m/s]

INTEGER, PARAMETER :: &
    hfDebug = 36  ! Diagnostic output

INTEGER :: &
    nerr  ! Error code

!------------------------------------------------------------------------------
!                               Namelists
!------------------------------------------------------------------------------
NAMELIST /testmmm_input_control/ &
    input_kind, npoints

NAMELIST /testmmm_input_1stkind/               &
    cmodel, cW20, cDBM, cETG, cMTM, cETGM,     &
    lprint, lW20, lDBM, lETG, lMTM, lETGM,     &
    rmin, rmaj, elong, ne, nh, nz, nf, zeff,   &
    te, ti, q, btor, zimp, aimp, ahyd, aimass, &
    wexbs, gne, gni, gnh, gnz, gte, gti, gq,   &
    gvtor, vtor, gvpol, vpol, gvpar, vpar     

!------------------------------------------------------------------------------
!                               Program Execution
!------------------------------------------------------------------------------
CALL SYSTEM_CLOCK(count_rate=count_rate, count_max=count_max)

! Open output.csv first, so that it's empty if errors occur before calling mmm
OPEN(hfOut, file='output.csv', form='formatted', status='replace', iostat=nerr)
IF (nerr /= 0) THEN
    PRINT '(A)', "ERROR: output.csv could not be opened for writing"
    STOP
END IF

OPEN(hfIn, file='input',  form='formatted', status='old', iostat=nerr)
IF (nerr /= 0) THEN
    PRINT '(A)', "ERROR: input file could not be opened for reading"
    STOP
END IF

READ(hfIn, NML=testmmm_input_control)
IF (npoints == BADINT) THEN
    PRINT '(A)', "ERROR: npoints for the number of radial points needs to be set"
    STOP
ELSE IF (input_kind /= 1) THEN
    PRINT '(A)', "ERROR: Unsupported input kind; please use testmmm"
    STOP
END IF

CALL initialize_arrays(npoints)
PRINT '(A)', "Input of the first kind (values) is detected. Processing..."
READ(hfIn, NML=testmmm_input_1stkind)
CLOSE(hfIn)

! Fill parameter arrays with default values from modmmm
CALL set_mmm_switches(cmmm=cswitch, lmmm=lswitch)

! Assign user specified parameters
DO i = 1, MAXNOPT 
    If (abs(cW20(i) - BADREAL) > 1E-6_R8) cswitch(i, KW20) = cW20(i)
    If (abs(cDBM(i) - BADREAL) > 1E-6_R8) cswitch(i, KDBM) = cDBM(i)
    If (abs(cMTM(i) - BADREAL) > 1E-6_R8) cswitch(i, KMTM) = cMTM(i)
    If (abs(cETG(i) - BADREAL) > 1E-6_R8) cswitch(i, KETG) = cETG(i)
    If (abs(cETGM(i) - BADREAL) > 1E-6_R8) cswitch(i, KETGM) = cETGM(i)

    if (lW20(i) /= BADINT) lswitch(i, KW20) = lW20(i)
    if (lDBM(i) /= BADINT) lswitch(i, KDBM) = lDBM(i)
    if (lMTM(i) /= BADINT) lswitch(i, KMTM) = lMTM(i)
    if (lETG(i) /= BADINT) lswitch(i, KETG) = lETG(i)
    if (lETGM(i) /= BADINT) lswitch(i, KETGM) = lETGM(i)
ENDDO

! Call and time mmm
CALL SYSTEM_CLOCK(tic)
CALL mmm(rmin=rmin, rmaj=rmaj, rmaj0=rmaj(1), elong=elong, ne=ne,       &
         nh=nh, nz=nz, nf=nf, zeff=zeff, te=te, ti=ti, q=q, btor=btor,  &
         zimp=zimp, aimp=aimp, ahyd=ahyd, aimass=aimass, wexbs=wexbs,   &
         gne=gne, gni=gni, gnh=gnh, gnz=gnz, gte=gte, gti=gti, gq=gq,   &
         gvtor=gvtor, vtor=vtor, gvpar=gvpar, vpol=vpol, gvpol=gvpol,   &       
         vpar=vpar, xti=xti, xdi=xdi, xte=xte, xdz=xdz, xvt=xvt,        &  
         xvp=xvp, xtiW20=xtiW20, xdiW20=xdiW20, xteW20=xteW20,          &
         xtiDBM=xtiDBM, xdiDBM=xdiDBM, xteDBM=xteDBM, xteETG=xteETG,    &   
         xteMTM=xteMTM, xdiETGM=xdiETGM, xteETGM=xteETGM, nerr=nerr,    &  
         gammaW20=gammaW20, omegaW20=omegaW20, gammaDBM=gammaDBM,       &
         omegaDBM=omegaDBM, gammaMTM=gammaMTM, omegaMTM=omegaMTM,       &                              
         gammaETGM=gammaETGM, omegaETGM=omegaETGM, dbsqprf=dbsqprf,     &
         nprout=hfDebug, cmodel=cmodel, npoints=npoints, lprint=lprint, &
         cswitch=cswitch, lswitch=lswitch, vconv=vconv, vflux=vflux)   
CALL SYSTEM_CLOCK(toc)

IF (nerr /= 0) THEN
    PRINT '(A, I3)', "ERROR: MMM finished with error code ", nerr
    STOP
END IF    

PRINT '(A, F13.6, A)', "MMM 8.2 finished successfully!  Run Time:", (toc - tic) / REAL(count_rate), "s"

! Write output variable names 
WRITE(hfOut,'("#"A11, 34A12)') &
    "rmin,", &
    "xti,", &
    "xdi,", &
    "xte,", &
    "xdz,", &
    "xvt,", &
    "xvp,", &
    "xtiW20,", &
    "xdiW20,", &
    "xteW20,", &
    "xtiDBM,", &
    "xdiDBM,", &
    "xteDBM,", &
    "xteETG,", &
    "xteMTM,", &
    "xteETGM,", &
    "xdiETGM,", &
    "gmaW20ii,", &
    "omgW20ii,", &
    "gmaW20ie,", &
    "omgW20ie,", &
    "gmaW20ei,", &
    "omgW20ei,", &
    "gmaW20ee,", &
    "omgW20ee,", &
    "gmaDBM,", &
    "omgDBM,", &
    "gmaMTM,", &
    "omgMTM,", &
    "gmaETGM,", &
    "omgETGM,", &
    "dbsqprf "

! Write output variable units
WRITE(hfOut,'("#"A11, 34 A12)') &
    "m,",                       &  ! rmin
    (/("m^2/s,", i=1, 16)/),    &  ! all diffusivities
    (/("s^-1,", i=1, 14)/),     &  ! all growth rates, frequencies
    ""                             ! dbsqprf

! Write output variable values
DO j = 1, npoints
    WRITE(hfOut,'(0P F11.6, A, 31(ES11.3, A))') &
        rmin(j), ',', &
        xti(j), ',', &
        xdi(j), ',', &
        xte(j), ',', & 
        xdz(j), ',', &
        xvt(j), ',', &
        xvp(j), ',', & 
        xtiW20(j), ',', &
        xdiW20(j), ',', &
        xteW20(j), ',', &
        xtiDBM(j), ',', &
        xdiDBM(j), ',', &
        xteDBM(j), ',', &
        xteETG(j), ',', &
        xteMTM(j), ',', &
        xteETGM(j), ',', &
        xdiETGM(j), ',', &
        gammaW20(1,j), ',', &  ! gmaW20ii 
        omegaW20(1,j), ',', &  ! omgW20ii 
        gammaW20(2,j), ',', &  ! gmaW20ie 
        omegaW20(2,j), ',', &  ! omgW20ie 
        gammaW20(3,j), ',', &  ! gmaW20ei 
        omegaW20(3,j), ',', &  ! omgW20ei 
        gammaW20(4,j), ',', &  ! gmaW20ee 
        omegaW20(4,j), ',', &  ! omgW20ee 
        gammaDBM(j), ',', &
        omegaDBM(j), ',', &
        gammaMTM(j), ',', &
        omegaMTM(j), ',', &
        gammaETGM(j), ',', &
        omegaETGM(j), ',', &
        dbsqprf(j)
END DO

CLOSE(hfOut)

!------------------------------------------------------------------------------
!                               Contains
!------------------------------------------------------------------------------

CONTAINS 
SUBROUTINE initialize_arrays(np)
    ! Note: Deallocation occurs naturally when the program ends

    INTEGER, INTENT(IN) :: np  ! Array dimension

    ! Input variables
    ALLOCATE(rmin(np)); rmin = BADREAL
    ALLOCATE(rmaj(np)); rmaj = BADREAL
    ALLOCATE(elong(np)); elong = BADREAL
    ALLOCATE(ne(np)); ne = BADREAL
    ALLOCATE(nh(np)); nh = BADREAL
    ALLOCATE(nz(np)); nz = BADREAL
    ALLOCATE(nf(np)); nf = BADREAL
    ALLOCATE(zeff(np)); zeff = BADREAL
    ALLOCATE(te(np)); te = BADREAL
    ALLOCATE(ti(np)); ti = BADREAL
    ALLOCATE(q(np)); q = BADREAL
    ALLOCATE(btor(np)); btor = BADREAL
    ALLOCATE(zimp(np)); zimp = BADREAL
    ALLOCATE(aimp(np)); aimp = BADREAL
    ALLOCATE(ahyd(np)); ahyd = BADREAL
    ALLOCATE(aimass(np)); aimass = BADREAL
    ALLOCATE(wexbs(np)); wexbs = BADREAL
    ALLOCATE(gne(np)); gne = BADREAL
    ALLOCATE(gni(np)); gni = BADREAL
    ALLOCATE(gnh(np)); gnh = BADREAL
    ALLOCATE(gnz(np)); gnz = BADREAL
    ALLOCATE(gte(np)); gte = BADREAL
    ALLOCATE(gti(np)); gti = BADREAL
    ALLOCATE(gq(np)); gq = BADREAL
    ALLOCATE(gvtor(np)); gvtor = BADREAL
    ALLOCATE(vtor(np)); vtor = BADREAL
    ALLOCATE(gvpol(np)); gvpol = BADREAL
    ALLOCATE(vpol(np)); vpol = BADREAL
    ALLOCATE(gvpar(np)); gvpar = BADREAL
    ALLOCATE(vpar(np)); vpar = BADREAL

    ! Output variables
    ALLOCATE(xti(np)); xti = 0_R8
    ALLOCATE(xdi(np)); xdi = 0_R8
    ALLOCATE(xte(np)); xte = 0_R8
    ALLOCATE(xdz(np)); xdz = 0_R8
    ALLOCATE(xvt(np)); xvt = 0_R8
    ALLOCATE(xvp(np)); xvp = 0_R8
    ALLOCATE(gammaDBM(np)); gammaDBM = 0_R8
    ALLOCATE(omegaDBM(np)); omegaDBM = 0_R8
    ALLOCATE(xtiW20(np)); xtiW20 = 0_R8
    ALLOCATE(xdiW20(np)); xdiW20 = 0_R8
    ALLOCATE(xteW20(np)); xteW20 = 0_R8
    ALLOCATE(xtiDBM(np)); xtiDBM = 0_R8
    ALLOCATE(xdiDBM(np)); xdiDBM = 0_R8
    ALLOCATE(xteDBM(np)); xteDBM = 0_R8
    ALLOCATE(xteETG(np)); xteETG = 0_R8
    ALLOCATE(xteETGM(np)); xteETGM = 0_R8
    ALLOCATE(xdiETGM(np)); xdiETGM = 0_R8
    ALLOCATE(gammaW20(4, np)); gammaW20 = 0_R8
    ALLOCATE(omegaW20(4, np)); omegaW20 = 0_R8
    ALLOCATE(gammaMTM(np)); gammaMTM = 0_R8
    ALLOCATE(omegaMTM(np)); omegaMTM = 0_R8
    ALLOCATE(gammaETGM(np)); gammaETGM = 0_R8
    ALLOCATE(omegaETGM(np)); omegaETGM = 0_R8
    ALLOCATE(dbsqprf(np)); dbsqprf = 0_R8
    ALLOCATE(xteMTM(np)); xteMTM = 0_R8
    ALLOCATE(vconv(6, np)); vconv = 0_R8
    ALLOCATE(vflux(6, np)); vflux = 0_R8
END SUBROUTINE initialize_arrays

END PROGRAM mmm_wrapper
