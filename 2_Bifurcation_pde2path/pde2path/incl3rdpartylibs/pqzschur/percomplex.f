#include "fintrf.h"
C  mex -v -g GCC='/usr/bin/gcc-4.4' percomplex.f -lmwlapack -lmwblas
C 
C PERCOMPLEX.F
C   Gateway function for complex periodic eigenvalue problems.
C
C   Periodic Hessenberg decomposition
C   task = 1 :        [Ar] = percomplex(1,A,s)
C                  [Qr,Ar] = percomplex(1,A,s)
C
C   Periodic Schur decomposition
C   task = 2 :        [Ar] = percomplex(2,A,s)
C                  [Qr,Ar] = percomplex(2,A,s)
C
C   Reordering of periodic Schur form:
C   task = 3 :        [Ar] = percomplex(3,A,s,select)
C                  [Qr,Ar] = percomplex(3,A,s,select)
C                  [Qr,Ar] = percomplex(3,Q,A,s,select)
C
C Purpose:
C   To solve complex periodic eigenvalue problems of the form
C
C                               S(2)                 S(K)
C          A(:,:,1)  *  A(:,:,2)     * ... * A(:,:,K),           (1)
C
C   where A is N-by-N-by-K and S is the signature array with values
C   1 or -1. 
C
C   task =  1: Computes the Hessenberg-triangular form, i.e., the matrix
C              A(:,:,1) is reduced to upper Hessenberg form while the
C              other matrices are triangularized.
C
C   task =  2: Computes the periodic Schur form, i.e., all matrices are
C              triangularized.
C
C   task =  3: Reorders the eigenvalues selected by the array SELECT
C              to the top left part of the generalized product (1).
C
C Input parameters: 
C   task   - integer option to determine the computation to perform as
C            described above.
C   A      - complex N-by-N*K matrix. A(:,(I-1)*N+1:I*N) corresponds to
C            the factor A(:,:,I) of (1). If task = 3 each A(:,:,I) must
C            be upper triangular.
C   Q      - complex N-by-N*K matrix, containing K unitary factors.
C   select - if task = 3: logical N-vector specifying the eigenvalues
C            to be reordered to the top part of the generalized product
C            (1).
C
C Output parameters:
C   Ar     - complex N-by-N*K matrix containing the updated factors of
C            (1).
C   Qr     - complex N-by-N*K matrix, containing K unitary factors.
C
C Contributor:
C   D. Kressner, Dresden, December 2003.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
      IMPLICIT NONE
C     .. Mex-file interface parameters ..
      mwPointer         PLHS(*), PRHS(*)
      INTEGER           NLHS, NRHS
    

C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Mex-file integer functions ..
      mwPointer           mxCreateFull, mxGetPi, mxGetPr
      mwPointer           mxGetM, mxGetN, mxIsComplex 
      integer mxIsNumeric
C
C     .. Scalar parameters used by subroutines ..
      INTEGER           INFO, K, LDA1, LDA2, LDQ1, LDQ2, LDWORK,
     $                  LZWORK, N
C
C     .. Allocatable arrays ..
C     !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: SELECT(:)
      INTEGER,          ALLOCATABLE :: S(:), SCAL(:)
      COMPLEX*16,       ALLOCATABLE :: A(:,:,:), ALPHA(:), BETA(:),
     $                                 Q(:,:,:), ZWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DWORK(:)
C     .. Local variables and constant dimension arrays ..
      LOGICAL           COMPQ, UPDQ
      CHARACTER         JOB
      CHARACTER*120     TEXT
      INTEGER           I, IP, M, NLM, NN, TASK
      DOUBLE PRECISION  TEMP 
C
C     .. External Subroutines ..
      EXTERNAL          ZLASET, ZPGEQZ, ZPGHRD, ZPGORD
C
C     ..Intrinsic Functions..
      INTRINSIC         DBLE, MAX, MIN
C
C     Check for proper number of arguments.
C
      
      IF ( NRHS.LT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'PERCOMPLEX requires at least 3 input arguments' )
      ELSE IF ( NLHS.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'PERCOMPLEX requires at most 2 output arguments' )
      END IF
C
C     Check dimensions of input parameters and read/set scalar
C     parameters.
C
      
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 )
     $   CALL mexErrMsgTxt( 'TASK must be a scalar' )
      
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR. 
     $     mxIsComplex( PRHS(1) ).EQ.1 )
     $   CALL mexErrMsgTxt( 'TASK must be an integer scalar' )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
     
      TASK = TEMP
C
      
      IF ( TASK.LT.1 .OR. TASK.GT.3 )
     $   CALL mexErrMsgTxt
     $        ( 'The only admissible values of TASK are 1, 2 or 3.')
C
      COMPQ  = .FALSE.
      UPDQ   = .FALSE.
C
      IF ( TASK.LE.2 ) THEN
         IF ( NLHS.GT.1 )
     $      COMPQ = .TRUE.
	   NLM = 2
	ELSE IF ( TASK.EQ.3 ) THEN
	   IF ( NLHS.GT.1 )
     $      COMPQ = .TRUE.
	   IF ( NRHS.GT.4 )
     $      UPDQ = .TRUE.
	   NLM = 2
	END IF
C
      IF ( NLHS.GT.NLM ) THEN
         WRITE( TEXT, '('' PERCOMPLEX requires at most '',I4,
     $            '' output arguments'')' ) NLM
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      IF ( ( TASK.LE.2 ) .AND. NRHS.GT.3 ) THEN
         CALL mexErrMsgTxt
     $           ( 'PERCOMPLEX requires at most 3 input arguments' )
      ELSE IF ( ( TASK.EQ.3 ) .AND. NRHS.GT.5 ) THEN
         CALL mexErrMsgTxt
     $           ( 'PERCOMPLEX requires at most 5 input arguments' )
	END IF 
C
      IF ( TASK.EQ.3 .AND. NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $           ( 'HAPACK_HAEIG requires at least 4 input arguments' )
	END IF
C
      IF ( UPDQ ) THEN
         IP = 3
      ELSE
         IP = 2
      END IF
C
C     Check A.
C
      N = mxGetM( PRHS(IP) )
      IF ( N.EQ.0 ) THEN
         K = 1
      ELSE
         K = mxGetN( PRHS(IP) )
         IF ( MOD( K, N ).NE.0 ) THEN
            WRITE( TEXT,
     $       '(''The number of columns of A must be a multiple of '',
     $       I7)' ) N
            CALL mexErrMsgTxt( TEXT )
         ELSE
            K = K / N
         END IF
      END IF
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR. 
     $     mxIsComplex( PRHS(IP) ).EQ.0 )
     $   CALL mexErrMsgTxt( 'A must be a complex matrix' )
	IP = IP + 1
C
C     Check S.
C
      M  = mxGetM( PRHS(IP) )
      NN = mxGetN( PRHS(IP) )
      IF ( .NOT.( ( M.EQ.1.AND.NN.EQ.K ).OR.
     $            ( M.EQ.K.AND.NN.EQ.1 ) ) ) THEN
         WRITE( TEXT, '(''S must be a vector of length  '', I7)' ) K
         CALL mexErrMsgTxt( TEXT )
	END IF
      IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR. 
     $     mxIsComplex( PRHS(IP) ).EQ.1 )
     $   CALL mexErrMsgTxt( 'S must be a real vector' )
      IP = IP + 1
C
C     Check SELECT.
C
      IF ( TASK.EQ.3 ) THEN
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( .NOT.( ( M.EQ.1.AND.NN.EQ.N ).OR.
     $               ( M.EQ.N.AND.NN.EQ.1 ) ) )
     $      CALL mexErrMsgTxt( 'SELECT must be a vector of length n' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR. 
     $        mxIsComplex( PRHS(IP) ).EQ.1 )
     $      CALL mexErrMsgTxt( 'SELECT must be a real vector' )
      END IF
C
      IF ( UPDQ ) THEN
C
C        Check Q.
C
         IP = 2
         M  = mxGetM( PRHS(IP) )
         NN = mxGetN( PRHS(IP) )
         IF ( M.NE.N .OR. NN.NE.K*N )
     $      CALL mexErrMsgTxt( 'Q must be an n-by-n-by-k matrix' )
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR. 
     $        mxIsComplex( PRHS(IP) ).EQ.0 )
     $   CALL mexErrMsgTxt( 'Q must be a complex matrix' )
      END IF
C
C     Determine the lenghts of working arrays.
C
      LDA1 = MAX( 1, N )
	LDA2 = MAX( 1, N )
	IF ( COMPQ .OR. UPDQ ) THEN
	   LDQ1 = MAX( 1, N )
	   LDQ2 = MAX( 1, N )
	ELSE
	   LDQ1 = 1
	   LDQ2 = 1
	END IF
C
C     Workspace requirements.
C
      IF ( TASK.EQ.1 ) THEN
	   LDWORK = MAX( 1, K, N )
	   LZWORK = MAX( 1, 2*N )
	ELSE IF ( TASK.EQ.2 ) THEN
	   LDWORK = MAX( 1, K, N )
	   LZWORK = MAX( 1, 2*N )
	ELSE IF ( TASK.EQ.3 ) THEN
	   LDWORK = MAX( 1, K, N )
	   LZWORK = MAX( 1, K, N )
	END IF
       

C
C     Allocate variable dimension local arrays.
C     !Fortran 90/95
C
	ALLOCATE ( A(LDA1,LDA2,K) )
	ALLOCATE ( S(K) )
	ALLOCATE ( DWORK(LDWORK) )
	ALLOCATE ( ZWORK(LZWORK) )
	IF ( COMPQ .OR. UPDQ )
     $   ALLOCATE ( Q(LDQ1,LDQ2,K) )
	IF ( TASK.GE.2 )
     $   ALLOCATE ( ALPHA(N), BETA(N), SCAL(N) )
	IF ( TASK.EQ.3 )
     $   ALLOCATE ( SELECT(N) )
C
C     Copy inputs from MATLAB workspace to locally allocated arrays.
C
      IP = 1
	IF ( UPDQ ) THEN
	   IP = IP + 1
	   CALL mxCopyPtrToComplex16( mxGetPr( PRHS(IP) ),
     $                           mxGetPi( PRHS(IP) ), Q, N*N*K )
	END IF
	IP = IP + 1
	CALL mxCopyPtrToComplex16( mxGetPr( PRHS(IP) ),
     $                           mxGetPi( PRHS(IP) ), A, N*N*K )
	IP = IP + 1
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DWORK, K )
	IF ( DWORK(1).NE.ONE )
     $   CALL mexErrMsgTxt( 'The first index, s(1), must be one' )
      DO 10  I = 1, K
	   IF ( ABS(DWORK(I)).NE.ONE )
     $      CALL mexErrMsgTxt( 'Each index, s(i), must be one or -one' )
         S(I) = INT( DWORK(I) )
   10 CONTINUE
      IF ( TASK.EQ.3 ) THEN
         IP = IP + 1
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), DWORK, N )
         DO 20  I = 1, N
            SELECT(I) = DWORK(I).EQ.ONE
   20    CONTINUE
      END IF
C
C     Do the actual computations.
C
      IF ( TASK.EQ.1 ) THEN
	   IF ( COMPQ ) THEN
	      JOB = 'I'
	   ELSE
	      JOB = 'N'
	   END IF
	   CALL ZPGHRD( JOB, K, N, 1, N, S, A, LDA1, LDA2, Q, LDQ1, LDQ2,
     $                DWORK, LDWORK, ZWORK, LZWORK, INFO )
      ELSE IF ( TASK.EQ.2 ) THEN
	   IF ( COMPQ ) THEN
	      JOB = 'I'
	   ELSE
	      JOB = 'N'
	   END IF
	   CALL ZPGHRD( JOB, K, N, 1, N, S, A, LDA1, LDA2, Q, LDQ1, LDQ2,
     $                DWORK, LDWORK, ZWORK, LZWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
	      IF ( COMPQ )
     $         JOB = 'V'
            CALL ZPGEQZ( 'Schur', JOB, K, N, 1, N, S, A, LDA1, LDA2,
     $                   ALPHA, BETA, SCAL, Q, LDQ1, LDQ2, DWORK,
     $                   LDWORK, ZWORK, LZWORK, INFO )
         END IF
	ELSE IF ( TASK.EQ.3 ) THEN
	   IF ( COMPQ .AND. .NOT.UPDQ ) THEN
	      DO 30 I = 1, K
	         CALL ZLASET( 'All', N, N, CZERO, CONE, Q(1,1,I), LDQ1 )
   30       CONTINUE
         END IF
	   CALL ZPGORD( COMPQ, K, N, S, SELECT, A, LDA1, LDA2,
     $                ALPHA, BETA, SCAL, Q, LDQ1, LDQ2, M,
     $                ZWORK, LZWORK, INFO )
      END IF
C
C     Copy output to MATLAB workspace.
C
      IP = 1
      IF ( INFO.EQ.0 ) THEN
	   IF ( COMPQ .AND. NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateFull( N, N*K, 1 )
            CALL mxCopyComplex16ToPtr( Q, mxGetPr( PLHS(IP) ),
     $                                 mxGetPi( PLHS(IP) ), N*N*K )
            IP = IP + 1
	   END IF
	   IF ( NLHS.GE.IP ) THEN
            PLHS(IP) = mxCreateFull( N, N*K, 1 )
            CALL mxCopyComplex16ToPtr( A, mxGetPr( PLHS(IP) ),
     $                                 mxGetPi( PLHS(IP) ), N*N*K )
            IP = IP + 1
	   END IF
      END IF
C
C     Deallocate local arrays.
C     !Fortran 90/95
C
	DEALLOCATE ( A, S, DWORK, ZWORK )
	IF ( COMPQ .OR. UPDQ )
     $   DEALLOCATE ( Q )
	IF ( TASK.GE.2 )
     $   DEALLOCATE ( ALPHA, BETA, SCAL )
	IF ( TASK.EQ.3 )
     $   DEALLOCATE ( SELECT )
C
C     Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         IF ( TASK.EQ.1 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,
     $             '' ON EXIT FROM ZPGHRD'')' )  INFO
         ELSE IF ( TASK.EQ.2 ) THEN
	      IF ( INFO.LT.0 ) THEN
               WRITE( TEXT, '('' INFO = '',I4,
     $                '' ON EXIT FROM ZPGHRD/ZPGEQZ'')' )  INFO
	      ELSE
               WRITE( TEXT,
     $         '('' The periodic QR algorithm failed to converge. '')' )
            END IF
         ELSE IF ( TASK.EQ.3 ) THEN
            WRITE( TEXT, '('' INFO = '',I4,
     $             '' ON EXIT FROM ZPGORD'')' )  INFO
	   END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of PERCOMPLEX ***
      END

! HUHU
      SUBROUTINE ZPGHRD( COMPQ, K, N, ILO, IHI, S, A, LDA1, LDA2, Q,
     $                   LDQ1, LDQ2, DWORK, LDWORK, ZWORK, LZWORK,
     $                   INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To reduce the general complex product
C
C                               S(2)                 S(K)
C          A(:,:,1)  *  A(:,:,2)     * ... * A(:,:,K)
C
C     to upper Hessenberg-triangular form, where A is N-by-N-by-K and S
C     is the signature array with values 1 or -1. The matrix A(:,:,1)
C     is reduced to upper Hessenberg form while the other matrices are
C     triangularized.
C     Simple and unblocked version.
C
C     If COMPQ = 'V' or COMPQ = 'I', then the unitary factors are
C     computed and stored in the array Q so that for S(I) = 1,
C
C                         H
C             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)
C                          H
C         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),
C
C     and for S(I) = -1,
C
C                                  H
C             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)
C                                   H
C         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).
C
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPQ   (input) CHARACTER*1
C             = 'N': do not modify Q.
C             = 'V': modify the array Q by the unitary transformations
C                    that are applied to the matrices in A to reduce them
C                    to Hessenberg-triangular form.
C             = 'I': like COMPQ='V', except that each matrix in Q will
C                    be initialized to the identity first.
C
C     Input/Output Parameters
C
C     K       (input) INTEGER
C             The number of matrices in A.  K >= 1.
C
C     N       (input) INTEGER
C             Order of each factor in A.  N >= 0.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that each factor in A is already upper
C             triangular in rows and columns 1:ILO-1 and IHI+1:N.
C             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
C
C     S       (input) INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array must
C             contain the factors of the general product to be reduced.
C             On exit, A(:,:,1) is overwritten by an upper Hessenberg
C             matrix and each A(:,:,I) for I not equal to 1 is
C             overwritten by an upper triangular matrix.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A. LDA1 >= max(1,N)
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A. LDA2 >= max(1,N)
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             If COMPQ='N': Q is not referenced.
C             If COMPQ='I': On entry, Q need not to be set, and on exit
C                           it contains the unitary transformations.
C             If COMPQ='V': On entry, Q must contain unitary matrices,
C                           and on exit this is overwritten by the
C                           updated transformations.
C
C     LDQ1    (input) INTEGER
C             The first leading dimension of Q. LDQ1 >= max(1,N)
C
C     LDQ2    (input) INTEGER
C             The second leading dimension of Q. LDQ2 >= max(1,N)
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the minimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N). 
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.  LZWORK >= MAX(1,2*N). 
C             For optimal performance this value should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C
C     METHOD
C
C     A slightly modified version of the periodic Hessenberg reduction
C     presented in [1] is used. For more details see [2].
C
C     REFERENCES
C
C     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P.
C         The periodic Schur decomposition; algorithm and applications.
C         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992.
C
C     [2] Kressner, D.
C         An efficient and reliable implementation of the periodic QZ
C         algorithm. In IFAC Workshop on Periodic Control Systems, 2001.
C
C     NUMERICAL ASPECTS
C
C     The implemented method is numerically backward stable.
C                                 3
C     The algorithm requires 0(K N ) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      CHARACTER*1       COMPQ
      INTEGER           K, ILO, IHI, INFO, LDA1, LDA2, LDQ1, LDQ2,
     $                  LDWORK, LZWORK, N
C     .. Array Arguments ..
      INTEGER           S(*)
      DOUBLE PRECISION  DWORK(*)
      COMPLEX*16        A(LDA1, LDA2, *), Q(LDQ1, LDQ2, *),
     $                  ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           INITQ, SOK, WANTQ
      INTEGER           ICOLS, IERR, IROWS, JCOL, JROW, L, WRKOPT
      DOUBLE PRECISION  CS
      COMPLEX*16        SN, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          ZGEQRF, ZGERQF, ZLACPY, ZLARTG, ZLASET, ZROT,
     $                  ZUNGQR, ZUNMQR, ZUNMRQ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCMPLX, DCONJG, INT, MAX
C
C     .. Executable Statements ..
C
      INFO   = 0
      WANTQ  = ( LSAME( COMPQ, 'V' ) ) .OR. ( LSAME( COMPQ, 'I' ) )
      INITQ  = LSAME( COMPQ, 'I' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.WANTQ .AND. .NOT.INITQ .AND. .NOT.LSAME(COMPQ,'N') ) THEN
         INFO = -1
      ELSE IF ( K .LT. 1 ) THEN
         INFO = -2
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE
         SOK = S(1).EQ.1
         DO 10  L = 2, K
            SOK = S(L).EQ.1 .OR. S(L).EQ.-1
   10    CONTINUE
         IF ( .NOT.SOK ) THEN
            INFO = -6
         ELSE IF ( LDA1 .LT. MAX(1, N) ) THEN
            INFO = -8
         ELSE IF ( LDA2 .LT. MAX(1, N) ) THEN
            INFO = -9
         ELSE IF ( WANTQ .AND. LDQ1.LT.MAX(1, N) ) THEN
            INFO = -11
         ELSE IF ( WANTQ .AND. LDQ2.LT.MAX(1, N) ) THEN
            INFO = -12
         ELSE IF ( LDWORK.LT.MAX(1,N) ) THEN
            INFO = -14
         ELSE IF ( LZWORK.LT.MAX(1,2*N) ) THEN
            INFO = -16
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'ZPGHRD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = DBLE(1)
         ZWORK(1) = CONE
         RETURN
      END IF
      WRKOPT = 2*N
C
C     Initialize Q if desired.
C
      IF ( INITQ )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ1 )
C
C     Transform A(2,:,:),...,A(K,:,:) to upper triangular form.
C
      DO 30  L = K, 2, -1
         IF ( S(L).EQ.1 ) THEN
C
C           Compute a QR Decomposition of A(:,:,L).
C
            IROWS = IHI + 1 - ILO
            ICOLS = N + 1 - ILO
            CALL ZGEQRF( IROWS, ICOLS, A(ILO,ILO,L), LDA1, ZWORK,
     $                   ZWORK(N+1), LZWORK-N, IERR )
            WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
C
C           Apply transformation to A(:,:,L-1).
C
            IF ( S(L-1).EQ.1 ) THEN
               CALL ZUNMQR( 'Right', 'No transpose', IHI, IROWS, IROWS,
     $                      A(ILO,ILO,L), LDA1, ZWORK, A(1,ILO,L-1),
     $                      LDA1, ZWORK(N+1), LZWORK-N, IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            ELSE
               CALL ZUNMQR( 'Left', 'Complex transpose', IROWS, ICOLS,
     $                      IROWS, A(ILO,ILO,L), LDA1, ZWORK, 
     $                      A(ILO,ILO,L-1), LDA1, ZWORK(N+1), LZWORK-N,
     $                      IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            END IF
C
C           Update transformation matrix Q(:,:,L).
C
            IF ( INITQ ) THEN
C               CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 )
               IF ( IROWS.GT.1 ) THEN
                  CALL ZLACPY( 'Lower', IROWS-1, IROWS-1,
     $                         A(ILO+1,ILO,L), LDA1, Q(ILO+1,ILO,L),
     $                         LDQ1 )
                  CALL ZUNGQR( IROWS, IROWS, IROWS, Q(ILO,ILO,L), LDQ1,
     $                         ZWORK, ZWORK(N+1), LZWORK-N, IERR)
                  WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
               END IF
            ELSE IF ( WANTQ ) THEN
               CALL ZUNMQR( 'Right', 'No transpose', N, IROWS, IROWS,
     $                       A(ILO,ILO,L), LDA1, ZWORK, Q(1,ILO,L),
     $                       LDQ1, ZWORK(N+1), LZWORK-N, IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            END IF
            IF ( IROWS.GT.1 )
     $         CALL ZLASET( 'Low', IROWS-1, IROWS-1, CZERO, CZERO,
     $                      A(ILO+1,ILO,L), LDA1)
         ELSE
C
C           Compute an RQ Decomposition of A(:,:,L).
C
            ICOLS = IHI + 1 - ILO
            CALL ZGERQF( IHI, ICOLS, A(1,ILO,L), LDA1, ZWORK,
     $                   ZWORK(N+1), LZWORK-N, IERR )
            WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
C
C           Apply transformation to A(:,:,L-1).
C
            IF ( S(L-1).EQ.1 ) THEN
               CALL ZUNMRQ( 'Right', 'Complex Transpose', IHI, ICOLS,
     $                      ICOLS, A(ILO,ILO,L), LDA1, ZWORK,
     $                      A(1,ILO,L-1), LDA1, ZWORK(N+1), LZWORK-N,
     $                      IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            ELSE
               CALL ZUNMRQ( 'Left', 'No transpose', ICOLS, N+1-ILO,
     $                      ICOLS, A(ILO,ILO,L), LDA1, ZWORK,
     $                      A(ILO,ILO,L-1), LDA1, ZWORK(N+1), LZWORK-N,
     $                      IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            END IF
C
C           Update transformation matrix Q(:,:,L).
C
            IF ( INITQ )
     $         CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 )
            IF ( INITQ.OR.WANTQ )
     $         CALL ZUNMRQ( 'Right', 'Complex transpose', N, ICOLS,
     $                      ICOLS, A(1,1,L), LDA1, ZWORK, Q(1,ILO,L),
     $                      LDQ1, ZWORK(N+1), LZWORK-N, IERR )
            WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            IF ( ICOLS.GT.1 )
     $         CALL ZLASET( 'Low', ICOLS-1, ICOLS-1, CZERO, CZERO,
     $                      A(ILO+1,ILO,L), LDA1 )
         END IF
   30 CONTINUE
C
C     Reduce A(:,:,1) to upper Hessenberg form.
C
      DO 110  JCOL = ILO, IHI - 2
C
C        Annihilate all elements below A(JCOL+1,JCOL,1).
C

         DO 40  JROW = IHI, JCOL + 2, -1
C
            TEMP = A(JROW-1,JCOL,1)
            CALL ZLARTG( TEMP, A(JROW,JCOL,1), CS, SN,
     $                   A(JROW-1,JCOL,1) )
            A(JROW,JCOL,1) = CZERO
            CALL ZROT( N-JCOL, A(JROW-1,JCOL+1,1), LDA1,
     $                 A(JROW,JCOL+1,1), LDA1, CS, SN )
            DWORK(JROW) = CS
            ZWORK(JROW) = SN
   40    CONTINUE
C
         IF ( WANTQ ) THEN
            DO 50  JROW = IHI, JCOL + 2, -1
               CALL ZROT( N, Q(1,JROW-1,1), 1, Q(1,JROW,1), 1,
     $                    DWORK(JROW), DCONJG( ZWORK(JROW) ) )
   50       CONTINUE
         END IF
C
C        Propagate transformations through A(:,:,K),...,A(:,:,2).
C
         DO 90  L = K, 2, -1
C
            IF ( S(L).EQ.1 ) THEN
               DO 60  JROW = IHI, JCOL + 2, -1
                  CALL ZROT( JROW, A(1,JROW-1,L), 1, A(1,JROW,L), 1,
     $                       DWORK(JROW), DCONJG( ZWORK(JROW) ) )
                  TEMP = A(JROW-1,JROW-1,L)
                  CALL ZLARTG( TEMP, A(JROW,JROW-1,L), CS, SN,
     $                         A(JROW-1,JROW-1,L) )
                  A(JROW,JROW-1,L) = CZERO
                  CALL ZROT( N-JROW+1, A(JROW-1,JROW,L), LDA1,
     $                       A(JROW,JROW,L), LDA1, CS, SN )
                  DWORK(JROW) = CS
                  ZWORK(JROW) = SN
   60          CONTINUE
            ELSE
               DO 70  JROW = IHI, JCOL + 2, -1
                  CALL ZROT( N+2-JROW, A(JROW-1,JROW-1,L), LDA1,
     $                       A(JROW,JROW-1,L), LDA1, DWORK(JROW),
     $                       ZWORK(JROW) )
                  TEMP = A(JROW,JROW,L)
                  CALL ZLARTG( TEMP, A(JROW,JROW-1,L), CS, SN,
     $                         A(JROW,JROW,L) )
                  A(JROW,JROW-1,L) = CZERO
                  CALL ZROT( JROW-1, A(1,JROW,L), 1, A(1,JROW-1,L), 1,
     $                       CS, SN )
                  DWORK(JROW) = CS
                  ZWORK(JROW) = -SN
   70          CONTINUE
            END IF
C
            IF ( WANTQ ) THEN
               DO 80  JROW = IHI, JCOL + 2, -1
                  CALL ZROT( N, Q(1,JROW-1,L), 1, Q(1,JROW,L), 1,
     $                       DWORK(JROW), DCONJG(ZWORK(JROW)) )
   80          CONTINUE
            END IF
   90    CONTINUE
C
C        Apply transformations to A(:,:,1).
C
         DO 100  JROW = IHI, JCOL + 2, -1
            CALL ZROT( IHI, A(1,JROW-1,L), 1, A(1,JROW,L), 1,
     $                 DWORK(JROW), DCONJG( ZWORK(JROW) ) )
  100    CONTINUE
  110 CONTINUE
      DWORK(1) = DBLE( N )
      ZWORK(1) = DCMPLX( WRKOPT, 0 )
      RETURN
C *** Last line of ZPGHRD ***
      END
	SUBROUTINE ZPGEX2( WANTQ, K, N, J, S, A, LDA1, LDA2, Q, LDQ1,
     $                   LDQ2, ZWORK, INFO )
	IMPLICIT NONE
C
C     PURPOSE
C
C     ZPGEX2 swaps adjacent diagonal 1-by-1 blocks in a complex
C     generalized matrix product,
C
C                                 S(2)                 S(K)
C          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K),
C
C     by unitary equivalence transformations. A must be in periodic
C     Schur form, that is, all factors of A must be upper triangular.
C
C     If WANTQ = .TRUE., then the unitary factors are computed and
C     stored in the array Q so that for S(I) = 1,
C
C                         H
C             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)
C                          H                                        (1)
C         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),
C
C     and for S(I) = -1,
C
C                                  H
C             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)
C                                   H                               (2)
C         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).
C
C     ARGUMEMTS
C
C     Mode Parameters
C
C     WANTQ   (input) LOGICAL
C             = .FALSE.: do not modify Q;
C             = .TRUE. : modify the array Q by the unitary
C                        transformations that are applied to the
C                        matrices in A for reordering.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of factors.  K >= 1.
C
C     N       (input)  INTEGER
C             The order of each factor in A.  N >= 2.
C
C     J       (input) INTEGER
C             The index of the first block to be swapped.  1 <= J < N.
C
C     S       (input)  INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array
C             must contain the factors in periodic Schur form
C             form, that is, all factors are upper triangular.
C             On exit, if INFO = 0, the leading N-by-N-by-K part of
C             this array contains the factors of the reordered periodic
C             Schur form.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A.  LDA1 >= MAX(1,N).
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A.  LDA2 >= MAX(1,N).
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             On entry, if WANTQ = .TRUE., the leading N-by-N-by-K part
C             of this array must contain the initial unitary factors
C             as described in (1)-(2).
C             On exit, if WANTQ = .TRUE., the leading N-by-N-by-K part
C             of this array contains the modified orthogonal factors as
C             described in (1)-(2).
C
C     LDQ1    (input)  INTEGER
C             The first leading dimension of Q.
C             If WANTQ = .TRUE.,  LDQ1 >= MAX(1,N).
C
C     LDQ2    (input)  INTEGER
C             The second leading dimension of Q.
C             If WANTQ = .TRUE.,  LDQ2 >= MAX(1,N).
C
C     Workspace
C
C     ZWORK   COMPLEX*16 array, dimension (K)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0       : succesful exit;
C             = 1       : the periodic QZ algorithm failed to converge.
C
C     METHOD
C
C     A complex version of the periodic QZ algorithm [1] with one
C     perfect shift is used. For more details see [2]. This is not
C     a safe method. It is advisable to check whether the eigenvalues
C     are really reorderd.
C
C     REFERENCES
C
C     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P.
C         The periodic Schur decomposition; algorithm and applications.
C         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992.
C
C     [2] Kressner, D.
C         An efficient and reliable implementation of the periodic QZ
C         algorithm. In IFAC Workshop on Periodic Control Systems, 2001.
C
C     NUMERICAL ASPECTS
C
C     The implemented method is numerically backward stable.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
	PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      LOGICAL           WANTQ
      INTEGER           J, K, INFO, LDA1, LDA2, LDQ1, LDQ2, N
C     .. Array Arguments ..
      INTEGER           S(*)
      COMPLEX*16        A(LDA1, LDA2, *), Q(LDQ1, LDQ2, *),
     $                  ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           USEZQ
      INTEGER           JITER, L, LN
      DOUBLE PRECISION  CS, CSF, RHS, SAFMAX, SAFMIN, SMLNUM, ULP
	COMPLEX*16        SN, SNF, TEMP
C     .. Local Arrays ..
      INTEGER           ISEED(4)
	COMPLEX*16        RND(3)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
	EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DLABAD, ZLARTG, ZROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, MAX
C
C     .. Executable Statements ..
C
C     Save perfect shift.
C
      DO 10  L = 1, K
	   ZWORK(L) = A(J,J,L)
   10 CONTINUE
      INFO = 0
      ISEED(1) = 1
	ISEED(2) = 0
	ISEED(3) = 0
	ISEED(4) = 0
      SAFMIN = DLAMCH( 'SafeMinimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Precision' )
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SAFMIN*( 2 / ULP )
C
C     If there are infinite or zero eigenvalues in critical positions,
C     then unshifted ZQ iterations must be used.
C
      USEZQ = .FALSE.
      DO 20  L = 1, K
         USEZQ = ( USEZQ.OR.
     $             ( ( S(L).EQ.1 ).AND.( A(J+1,J+1,L).EQ.CZERO ) ).OR.
     $             ( ( S(L).EQ.-1 ).AND.( A(J,J,L).EQ.CZERO ) ) )
   20 CONTINUE
C
C     Destroy any triangular structure.
C
      CALL ZLARNV( 2, ISEED, 2, RND )
      CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
      CSF = CS
	SNF = SN
      DO 30  L = 1, K
	   IF ( WANTQ ) THEN
            CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS, DCONJG( SN ) )
	   END IF
         IF ( S(L).EQ.1 ) THEN
	      CALL ZROT( N-J+1, A(J,J,L), LDA1, A(J+1,J,L), LDA1, CS, SN )
	   ELSE
	      CALL ZROT( J+1, A(1,J,L), 1, A(1,J+1,L), 1, CS,
     $                 DCONJG( SN ) )	      
	   END IF
	   IF ( L.EQ.K ) THEN
	      CS = CSF
	      SN = SNF
	   ELSE
            CALL ZLARNV( 2, ISEED, 2, RND )
            CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
	   END IF
	   IF ( S(L).EQ.1 ) THEN
	      CALL ZROT( J+1, A(1,J,L), 1, A(1,J+1,L), 1, CS,
     $                 DCONJG( SN ) )
	   ELSE
            CALL ZROT( N-J+1, A(J,J,L), LDA1, A(J+1,J,L), LDA1, CS, SN )
	   END IF
   30 CONTINUE
C
      IF ( USEZQ ) THEN
	   DO 50  JITER = 1, 10
            DO 40  L = 1, K
	         IF ( S(L).EQ.1 ) THEN
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J, A(1,J+1,L), 1, A(1,J,L), 1, CS, SN )
	            SN = -SN
	         ELSE
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( N-J, A(J,J+1,L), LDA1, A(J+1,J+1,L), LDA1,
     $                       CS, SN )
	         END IF
	         LN = L + 1
	         IF ( LN.GT.K )  LN = 1
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,LN), 1, Q(1,J+1,LN), 1, CS,
     $                       DCONJG( SN ) )
               END IF
	         IF ( S(LN).EQ.1 ) THEN
                  CALL ZROT( N-J+1, A(J,J,LN), LDA1, A(J+1,J,LN), LDA1,
     $                       CS, SN )
	         ELSE
                  CALL ZROT( J+1, A(1,J,LN), 1, A(1,J+1,LN), 1, CS,
     $                       DCONJG( SN ) )
	         END IF
   40       CONTINUE
            RHS = MAX( ABS( A(J,J,1) ), ABS( A(J+1,J+1,1) ) )
            IF ( RHS.EQ.ZERO )
     $         RHS = ABS( A(J,J+1,1) )
	      RHS = MAX( ULP*RHS, SMLNUM )
            IF ( ABS( A(J+1,J,1) ).LE.RHS ) THEN
	         A(J+1,J,1) = CZERO
               GO TO 90
            END IF
   50    CONTINUE
         INFO = 1
	   GO TO 90
	END IF

      DO 80  JITER = 1, 42
C
C        Complex single shift.
C
         IF ( ( JITER / 10 )*10.EQ.JITER ) THEN
C
C           Random/exceptional shift.
C
            CALL ZLARNV( 2, ISEED, 2, RND )
            CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
         ELSE
	      CALL ZLARTG( CONE, CONE, CS, SN, TEMP )
            DO 60  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
	            CALL ZLARTG( A(J,J,L)*CS, ZWORK(L)*DCONJG(SN),
     $                         CS, SN, TEMP )
               ELSE
	            CALL ZLARTG( ZWORK(L)*CS, -A(J,J,L)*DCONJG(SN),
     $                         CS, SN, TEMP )
	            SN = -SN
               END IF
   60       CONTINUE
            CALL ZLARTG( A(J,J,1)*CS - DCONJG(SN)*ZWORK(1),
     $                   A(J+1,J,1)*CS, CS, SN, TEMP )
         END IF
C
C        Do one QZ sweep.
C
         CALL ZROT( N-J+1, A(J,J,1), LDA1, A(J+1,J,1), LDA1,
     $              CS, SN )
         IF ( WANTQ ) THEN
            CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, CS,
     $                 DCONJG( SN ) )
         END IF
C
C        Propagate rotation through AK, ..., A2 to A1.
C
         DO 70  L = K, 2, -1
            IF ( S(L).EQ.1 ) THEN
               CALL ZROT( J+1, A(1,J,L), 1, A(1,J+1,L), 1, CS,
     $                    DCONJG( SN ) )
               TEMP = A(J,J,L)
               CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
               A(J+1,J,L) = CZERO
               CALL ZROT( N-J, A(J,J+1,L), LDA1, A(J+1,J+1,L), LDA1,
     $                    CS, SN )
            ELSE
               CALL ZROT( N-J+1, A(J,J,L), LDA1, A(J+1,J,L), LDA1, CS,
     $                    SN )
               TEMP = A(J+1,J+1,L)
               CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
               A(J+1,J,L) = CZERO
               CALL ZROT( J, A(1,J+1,L), 1, A(1,J,L), 1, CS, SN )
               SN = -SN
            END IF
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS,
     $                    DCONJG( SN ) )
            END IF
   70    CONTINUE
         CALL ZROT( J+1, A(1,J,1), 1, A(1,J+1,1), 1, CS, DCONJG( SN ) )
C
C        Test for deflation.
C
         RHS = MAX( ABS( A(J,J,1) ), ABS( A(J+1,J+1,1) ) )
         IF ( RHS.EQ.ZERO )
     $      RHS = ABS( A(J,J+1,1) )
	   RHS = MAX( ULP*RHS, SMLNUM )
         IF ( ABS( A(J+1,J,1) ).LE.RHS ) THEN
	      A(J+1,J,1) = CZERO
            GO TO 90
         END IF
   80 CONTINUE
C     
C     Not converged.
C
      INFO = 1
C
   90 CONTINUE
C
C     Check for singular triangular factors.
C
      DO 100  L = 1, K
         RHS = MAX( ABS( A(J,J+1,L) ), ABS( A(J+1,J+1,L) ) )
	   RHS = MAX( ULP*RHS, SMLNUM )
         IF ( ABS( A(J,J,L) ).LE.RHS ) THEN
	      A(J,J,L) = CZERO
         END IF
         RHS = MAX( ABS( A(J,J+1,L) ), ABS( A(J,J,L) ) )
	   RHS = MAX( ULP*RHS, SMLNUM )
         IF ( ABS( A(J+1,J+1,L) ).LE.RHS ) THEN
	      A(J+1,J+1,L) = CZERO
         END IF
  100 CONTINUE
      RETURN
C *** Last line of ZPGEX2 ***
      END
      SUBROUTINE ZLAPR1( BASE, K, S, A, INCA, ALPHA, BETA, SCAL )
      IMPLICIT NONE
C
C     PURPOSE
C
C     Computes the general product of K complex scalars trying to avoid
C     over- and underflow.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     BASE    (input)  DOUBLE PRECISION
C             Machine base.
C
C     K       (input)  INTEGER
C             The number of scalars.  K >= 0.
C
C     S       (input)  INTEGER array, dimension (K)
C             The signature array. Each entry of S must be 1 or -1.
C
C     A       (input)  COMPLEX*16 array, dimension (K)
C             Vector of real scalars.
C
C     INCA    (input)  INTEGER
C             Increment for the array A. incA <> 0.
C
C     ALPHA   (output)  COMPLEX*16
C             ALPHA is a real scalar with 1.0 <= ABS(ALPHA) < BASE such
C             that
C
C                  ALPHA / BETA * BASE**(SCAL)
C
C             is the general product of A. 
C
C     BETA    (output)  COMPLEX*16
C             BETA is either 0.0 or 1.0.
C             See also the description of ALPHA.
C
C     SCAL    (output)  INTEGER
C             Scaling factor, see ALPHA.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
	PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Array Arguments ..
      INTEGER           S(*)
      COMPLEX*16        A(*)
C     .. Scalar Arguments ..
      INTEGER           INCA, K, SCAL
	DOUBLE PRECISION  BASE
      COMPLEX*16        ALPHA, BETA
C     .. Local Scalars ..
      INTEGER           I, INDA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX
C
C     .. Executable Statements ..
C
      ALPHA = CONE
      BETA = CONE
      SCAL = 0
	INDA = 1
	DO 50  I = 1, K
	   IF ( S(I).EQ.1 ) THEN
	      ALPHA = ALPHA * A(INDA)
         ELSE
            IF ( A(INDA).EQ.CZERO ) THEN
	         BETA = CZERO
	      ELSE
	         ALPHA = ALPHA / A(INDA)
	      END IF
	   END IF
	   IF ( ABS( ALPHA ).EQ.ZERO ) THEN
	      SCAL = 0
	      ALPHA = CZERO
         ELSE
   10	      IF ( ABS( ALPHA ).GE.ONE )  GO TO 20
               ALPHA = DCMPLX( BASE, ZERO )*ALPHA
	         SCAL = SCAL - 1
	      GO TO 10
   20       CONTINUE
   30       IF ( ABS( ALPHA ).LT.BASE ) GO TO 40
               ALPHA = ALPHA / DCMPLX( BASE, ZERO )
	         SCAL = SCAL + 1
            GO TO 30
   40       CONTINUE
	   END IF
	   INDA = INDA + INCA
   50 CONTINUE
      RETURN
C *** Last line of ZLAPR1 ***
      END
      SUBROUTINE ZPGEQZ( JOB, COMPQ, K, N, ILO, IHI, S, A, LDA1, LDA2, 
     $                   ALPHA, BETA, SCAL, Q, LDQ1, LDQ2, DWORK,
     $                   LDWORK, ZWORK, LZWORK, INFO )
	IMPLICIT NONE
C
C     PURPOSE
C
C     ZPGEQZ implements a single-shift version of the periodic QZ
C     method for finding the eigenvalues of the complex generalized
C     matrix product
C
C                                 S(2)                 S(K)
C          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K).
C
C     In addition, A may be reduced to periodic Schur form by unitary
C     transformations: all factors A(:,:,i) become upper triangular.
C
C     If COMPQ = 'V' or COMPZ = 'I', then the unitary factors are
C     computed and stored in the array Q so that for S(I) = 1,
C
C                         H
C             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)
C                          H                                        (1)
C         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),
C
C     and for S(I) = -1,
C
C                                  H
C             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)
C                                   H                               (2)
C         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).
C
C     ARGUMEMTS
C
C     Mode Parameters
C
C     JOB     (input) CHARACTER*1
C             = 'E': compute only the eigenvalues; A will not
C                    necessarily be put into periodic Schur form.
C             = 'S': put A into periodic Schur form, as well
C                    as computing the eigenvalues contained in ALPHAR,
C                    ALPHAI, BETA and SCAL.
C
C     COMPQ   (input) CHARACTER*1
C             = 'N': do not modify Q.
C             = 'V': modify the array Q by the unitary transformations
C                    that are applied to the matrices in A to reduce them
C                    to periodic Schur form.
C             = 'I': like COMPQ='V', except that each matrix in Q will
C                    be initialized to the identity first.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of factors.  K >= 1.
C
C     N       (input)  INTEGER
C             The order of each factor in A.  N >= 0.
C
C     ILO     (input)  INTEGER
C     IHI     (input)  INTEGER
C             It is assumed that each factor in A is already upper
C             triangular in rows and columns 1:ILO-1 and IHI+1:N.
C             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
C
C     S       (input)  INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array
C             must contain the factors in upper Hessenberg-triangular
C             form, that is, A(:,:,1) is upper Hessenberg and the other
C             factors are upper triangular.
C             On exit, if JOB = 'S' and INFO = 0, the leading
C             N-by-N-by-K part of this array contains the factors of
C             A in periodic Schur form. All factors are reduched to
C             upper triangular form and, moreover, A(:,:,2),...,
C             A(:,:,K) are normalized so that their diagonals contain
C             nonnegative real numbers.
C             On exit, if JOB = 'E', then the leading N-by-N-by-K part
C             of this array contains meaningless elements.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A.  LDA1 >= MAX(1,N).
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A.  LDA2 >= MAX(1,N).
C
C     ALPHA   (output) COMPLEX*16 array, dimension (N)
C             On exit, if INFO = 0, the leading N elements of this
C             array contain the scaled eigenvalues of A. The i-th
C             eigenvalue of A is given by
C
C             ALPHA(I) / BETA(I) * BASE**(SCAL(I)),
C
C             where 1.0 <= ABS(ALPHA(I)) < BASE and BASE is the machine
C             base (normally 2.0).
C
C     BETA    (output) COMPLEX*16 array, dimension (N)
C             On exit, if INFO = 0, the leading N elements of this
C             array contain indicators for infinite eigenvalues. That
C             is, if BETA(I) = 0.0, then the i-th eigenvalue is
C             infinite. Otherwise BETA(I) is set to 1.0.
C
C     SCAL    (output) INTEGER array, dimension (N)
C             On exit, if INFO = 0, the leading N elements of this
C             array contain the scaling parameters for the eigenvalues
C             of A.
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             On entry, if COMPQ = 'V', the leading N-by-N-by-K part
C             of this array must contain the initial unitary factors
C             as described in (1)-(2).
C             On exit, if COMPQ = 'V' or COMPQ = 'I', the leading
C             N-by-N-by-K part of this array contains the modified
C             orthogonal factors as described in (1)-(2).
C
C     LDQ1    (input)  INTEGER
C             The first leading dimension of Q.  LDQ1 >= MAX(1,N).
C
C     LDQ2    (input)  INTEGER
C             The second leading dimension of Q.  LDQ2 >= MAX(1,N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the minimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N). 
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the minimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.  LZWORK >= MAX(1,N). 
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0       : succesful exit;
C             < 0       : if INFO = -i, the i-th argument had an
C                         illegal value;
C             = 1,..,N  : the periodic QZ iteration did not converge.
C                         A is not in periodic Schur form, but
C                         ALPHA(I), BETA(I) and SCAL(I), for
C                         I = INFO+1,...,N should be correct.
C
C     METHOD
C
C     A slightly modified version of the periodic QZ algorithm is
C     used. For more details see [2].
C
C     REFERENCES
C
C     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P.
C         The periodic Schur decomposition; algorithm and applications.
C         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992.
C
C     [2] Kressner, D.
C         An efficient and reliable implementation of the periodic QZ
C         algorithm. In IFAC Workshop on Periodic Control Systems, 2001.
C
C     NUMERICAL ASPECTS
C
C     The implemented method is numerically backward stable.
C                                 3
C     The algorithm requires 0(K N ) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      CHARACTER*1       COMPQ, JOB
      INTEGER           K, ILO, IHI, INFO, LDA1, LDA2, LDQ1, LDQ2,
     $                  LDWORK, LZWORK, N
C     .. Array Arguments ..
      INTEGER           S(*), SCAL(*)
	DOUBLE PRECISION  DWORK(*)
      COMPLEX*16        A(LDA1, LDA2, *), ALPHA(*), BETA(*),
     $                  Q(LDQ1, LDQ2, *), ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           INITQ, LSCHR, WANTQ, SOK
      INTEGER           IFIRST, IFRSTM, IITER, ILAST, ILASTM, IN, J, J1,
     $                  JDEF, JITER, JLO, L, LDEF, LN, MAXIT, NTRA,
     $                  ZITER
      DOUBLE PRECISION  ABST, BASE, CS, SAFMIN, SAFMAX, SMLNUM, ULP,
     $                  TOL
	COMPLEX*16        SN, TEMP
C     .. Local Arrays ..
      INTEGER           ISEED(4)
	COMPLEX*16        RND(4)
C     .. External Functions ..
      LOGICAL           LSAME
	DOUBLE PRECISION  DLAMCH, ZLANHS
      EXTERNAL          DLAMCH, LSAME, ZLANHS
C     .. External Subroutines ..
      EXTERNAL          DLABAD, ZLAPR1, ZLARNV, ZLARTG, ZROT, ZSCAL
      INTRINSIC         ABS, DBLE, DCMPLX, DCONJG, INT, LOG, MAX
C
C     .. Executable Statements ..
C
      INFO = 0
      LSCHR = LSAME( JOB,'S' )
      WANTQ = LSAME( COMPQ,'V' ).OR.LSAME( COMPQ,'I' )
      INITQ = LSAME( COMPQ,'I' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT. ( LSCHR .OR. LSAME( JOB,'E' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTQ .AND. .NOT.INITQ .AND.
     $         .NOT.LSAME(COMPQ,'N') ) THEN
         INFO = -2
      ELSE IF ( K.LT.1 ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF ( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN   
         INFO = -6
      ELSE
	   SOK = S(1).EQ.1
	   DO 10  L = 2, K
	      SOK = S(L).EQ.1 .OR. S(L).EQ.-1
   10    CONTINUE
         IF ( .NOT.SOK ) THEN
	      INFO = -7
         ELSE IF ( LDA1 .LT. MAX(1, N) ) THEN
            INFO = -9
         ELSE IF ( LDA2 .LT. MAX(1, N) ) THEN	
            INFO = -10
         ELSE IF ( WANTQ .AND. LDQ1.LT.MAX(1, N) ) THEN
            INFO = -15
         ELSE IF ( WANTQ .AND. LDQ2.LT.MAX(1, N) ) THEN
            INFO = -16
         ELSE IF ( LDWORK.LT.MAX(1,N) ) THEN
            INFO = -18
	   ELSE IF ( LZWORK.LT.MAX(1,N) ) THEN
            INFO = -20
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPGEQZ', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         ZWORK(1) = CONE
         RETURN
      END IF
C
C     Initialize Q.
C
      IF ( INITQ ) THEN
         DO 20  L = 1, K
            CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 )
   20    CONTINUE
      END IF
C
C     Machine Constants
C
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'SafeMinimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Precision' )
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SAFMIN*( IN / ULP )
      BASE = DLAMCH( 'Base' )
      IF ( K.GE.INT( LOG( DLAMCH('Underflow') ) / LOG( ULP ) ) ) THEN
C
C        Start Iteration with a controlled zero shift.
C
         ZITER = -1
      ELSE
         ZITER = 0
      END IF
C
C     Set Eigenvalues IHI+1:N
C
      DO 30  J = IHI + 1, N
	   CALL ZLAPR1( BASE, K, S, A(J,J,1), LDA1*LDA2, ALPHA(J),
     $                BETA(J), SCAL(J) )
   30 CONTINUE
C
C     If IHI < ILO, skip QZ steps
C
      If ( IHI.LT.ILO )  GOTO 470
C
C     MAIN PERIODIC QZ ITERATION LOOP
C     
C     Initialize dynamic indices
C
C     Eigenvalues ILAST+1:N have been found.
C        Column operations modify rows IFRSTM:whatever.
C        Row operations modify columns whatever:ILASTM.
C
C     If only eigenvalues are being computed, then
C        IFRSTM is the row of the last splitting row above row ILAST;
C        this is always at least ILO.
C     IITER counts iterations since the last eigenvalue was found,
C        to tell when to use an observed zero or random shift.
C     MAXIT is the maximum number of QZ sweeps allowed.
C
      ILAST = IHI
      IF ( LSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ISEED(1) = 1
      ISEED(2) = 0
      ISEED(3) = 0
      ISEED(4) = 0
      MAXIT = 30 * IN
C
      DO  460 JITER = 1, MAXIT
C
C        Special Case: ILAST = ILO
C
         IF ( ILAST.EQ.ILO )  GOTO 400
C
C        **************************************************************
C        *                     CHECK FOR DEFLATION                    *
C        **************************************************************
C
C        Test 1:  Deflation in the Hessenberg matrix.
C
         JLO = ILO
         DO 40  J = ILAST, ILO+1, -1
            TOL = ABS( A(J-1,J-1,1) ) + ABS( A(J,J,1) )
            IF ( TOL.EQ.ZERO )
     $         TOL = ZLANHS( '1',J-ILO+1, A(ILO,ILO,1), LDA1,
     $                       DWORK(1) )
            TOL = MAX( ULP*TOL, SMLNUM )
            IF ( ABS( A(J,J-1,1) ).LE.TOL ) THEN
               A(J,J-1,1) = CZERO
               JLO = J
               IF ( J.EQ.ILAST )  GOTO 400
               GOTO 50
            END IF
   40    CONTINUE
C
   50    CONTINUE
C
C        Test 2:  Deflation in the triangular matrices with index 1.
C
         DO 70  LDEF = 2,K
            IF ( S(LDEF).EQ.1 ) THEN
               DO 60  J = ILAST, JLO, -1
                  IF ( J.EQ.ILAST ) THEN
                     TOL = ABS( A(J-1,J,LDEF) )
                  ELSE IF ( J.EQ.JLO ) THEN
                     TOL = ABS( A(J,J+1,LDEF) )
                  ELSE
                     TOL = ABS( A(J-1,J,LDEF) ) + ABS( A(J,J+1,LDEF) )
                  END IF
                  IF ( TOL.EQ.ZERO )
     $               TOL = ZLANHS( '1', J-JLO+1, A(JLO,JLO,LDEF), LDA1,
     $                             DWORK(1) )
                  TOL = MAX( ULP*TOL, SMLNUM )
                  IF ( ABS( A(J,J,LDEF) ).LE.TOL ) THEN
                     A(J,J,LDEF) = CZERO
                     GOTO 180
                  END IF          
   60          CONTINUE
            END IF
   70    CONTINUE
C
C        Test 3:  Deflation in the triangular matrices with index -1.
C
         DO 90  LDEF = 2,K
            IF ( S(LDEF).EQ.-1 ) THEN
               DO 80  J = ILAST, JLO, -1
                  IF ( J.EQ.ILAST ) THEN
                     TOL = ABS( A(J-1,J,LDEF) )
                  ELSE IF ( J.EQ.JLO ) THEN
                     TOL = ABS( A(J,J+1,LDEF) )
                  ELSE
                     TOL = ABS( A(J-1,J,LDEF) ) + ABS( A(J,J+1,LDEF) )
                  END IF
                  IF ( TOL.EQ.ZERO )
     $                TOL = ZLANHS( '1', J-JLO+1, A(JLO,JLO,LDEF), LDA1,
     $                              DWORK(1) )
                  TOL = MAX( ULP*TOL, SMLNUM )
                  IF ( ABS( A(J,J,LDEF) ).LE.TOL ) THEN
                     A(J,J,LDEF) = CZERO
                     GOTO 330
                  END IF          
   80          CONTINUE
            END IF
   90    CONTINUE
C
C        Test 4:  Controlled zero shift.
C
         IF ( ZITER.GE.7 .OR. ZITER.LT.0 ) THEN
C
C           Make Hessenberg matrix upper triangular.
C
            DO 100 J = JLO, ILAST-1
               TEMP = A(J,J,1)
               CALL ZLARTG( TEMP, A(J+1,J,1), CS, SN, A(J,J,1) )
               A(J+1,J,1) = CZERO
               CALL ZROT( ILASTM-J, A(J,J+1,1), LDA1,
     $                    A(J+1,J+1,1), LDA1, CS, SN )
               DWORK(J) = CS
               ZWORK(J) = SN
  100       CONTINUE
            IF ( WANTQ ) THEN
               DO 110  J = JLO, ILAST-1
                  CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, DWORK(J),
     $                      DCONJG( ZWORK(J) ) )
  110          CONTINUE
            END IF
C
C           Propagate Transformations back to A_1.
C
            DO 150  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
                  DO 120  J = JLO, ILAST-1
                     CS = DWORK(J)
                     SN = ZWORK(J)
                     IF ( SN.NE.CZERO ) THEN
                        CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                             A(IFRSTM,J+1,L), 1, CS,
     $                             DCONJG( SN ) )
C
C                       Check for deflation
C
                        TOL = ABS( A(J,J,L) ) + ABS( A(J+1,J+1,L) )
                        IF ( TOL.EQ.ZERO )
     $                     TOL = ZLANHS( '1',J-JLO+2, A(JLO,JLO,L),
     $                                   LDA1, DWORK(1) )
                        TOL = MAX( ULP*TOL, SMLNUM )
                        IF ( ABS( A(J+1,J,L) ).LE.TOL ) THEN
                           CS = ONE
                           SN = CZERO
                           A(J+1,J,L) = CZERO
                        END IF
C
                        TEMP = A(J,J,L)
                        CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN,
     $                               A(J,J,L) )
                        A(J+1,J,L) = CZERO
                        CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                             A(J+1,J+1,L), LDA1, CS, SN )
                     END IF
                     DWORK(J) = CS
                     ZWORK(J) = SN
  120             CONTINUE
               ELSE
                  DO 130  J = JLO, ILAST-1
                     CS = DWORK(J)
                     SN = ZWORK(J)
                     IF ( SN.NE.CZERO ) THEN
                        CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1,
     $                             A(J+1,J,L), LDA1, CS, SN )
C
C                       Check for deflation
C
                        TOL = ABS( A(J,J,L) ) + ABS( A(J+1,J+1,L) )
                        IF ( TOL.EQ.ZERO )
     $                     TOL = ZLANHS( '1',J-JLO+2, A(JLO,JLO,L),
     $                                   LDA1, DWORK(1) )
                        TOL = MAX( ULP*TOL, SMLNUM )
                        IF ( ABS( A(J+1,J,L) ).LE.TOL ) THEN
                           CS = ONE
                           SN = CZERO
                           A(J+1,J,L) = CZERO
                        END IF
C
                        TEMP = A(J+1,J+1,L)
                        CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN,
     $                               A(J+1,J+1,L) )
                        A(J+1,J,L) = CZERO
                        CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                             A(IFRSTM,J,L), 1, CS, SN )
                        DWORK(J) = CS
                        ZWORK(J) = -SN
                     END IF
  130             CONTINUE
               END IF
C
               IF ( WANTQ ) THEN
                  DO 140  J = JLO, ILAST-1
                     CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, DWORK(J),
     $                          DCONJG( ZWORK(J) ) )
  140             CONTINUE
               END IF
  150       CONTINUE

C
C           Apply the transformations to the right hand side of the
C           Hessenberg factor.
C
            ZITER = 0
            DO 160  J = JLO, ILAST-1
               CS = DWORK(J)
               SN = ZWORK(J)
               CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,1), 1,
     $                    A(IFRSTM,J+1,1), 1, CS, DCONJG( SN ) )
	         IF ( SN.EQ.CZERO )
     $            ZITER = 1
  160       CONTINUE
C
C           No QZ iteration.
C
            GOTO 450
         END IF
C
C        **************************************************************
C        *                     HANDLE DEFLATIONS                      *
C        **************************************************************
C
C        Case I: Deflation occurs in the Hessenberg matrix. The QZ
C                iteration is only applied to the JLO:ILAST part.
C
  170    CONTINUE
         IFIRST = JLO
C
C        Go to the periodic QZ steps
C
         GOTO 410
C
C        Case II: Deflation occurs in a triangular matrix with index 1.
C
C        Do an unshifted periodic QZ step.
C
  180    JDEF = J
         DO 190  J = JLO, JDEF-1
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J+1,J,1), CS, SN, A(J,J,1) )
            A(J+1,J,1) = CZERO
            CALL ZROT( ILASTM-J, A(J,J+1,1), LDA1, A(J+1,J+1,1), LDA1,
     $                 CS, SN )
            DWORK(J) = CS
            ZWORK(J) = SN
  190    CONTINUE
         IF ( WANTQ ) THEN
            DO 200  J = JLO, JDEF-1
               CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, DWORK(J),
     $                    DCONJG( ZWORK(J) ) )
  200       CONTINUE
         END IF
C
C        Propagate the transformations through the triangular matrices.
C        Due to the zero element on the diagonal of the LDEF-th
C        factor the number of transformations drops by one.
C
         DO 240  L = K, 2, -1
            IF ( L.LT.LDEF ) THEN
               NTRA = JDEF-2
            ELSE
               NTRA = JDEF-1
            END IF
            IF ( S(L).EQ.1 ) THEN
               DO 210  J = JLO, NTRA
                  CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J+1,L), 1, DWORK(J),
     $                       DCONJG( ZWORK (J) ) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                       A(J+1,J+1,L), LDA1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = SN
  210          CONTINUE
            ELSE
               DO 220  J = JLO, NTRA
                  CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1, A(J+1,J,L),
     $                       LDA1, DWORK(J), ZWORK(J) )
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = -SN
  220          CONTINUE
            END IF
            IF ( WANTQ ) THEN
               DO 230  J = JLO, NTRA
                  CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, DWORK(J),
     $                       DCONJG( ZWORK(J) ) )
  230          CONTINUE
            END IF
  240    CONTINUE
C
C        Apply the transformations to the right hand side of the
C        Hessenberg factor.
C
         DO 250  J = JLO, JDEF-2
            CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,1), 1, A(IFRSTM,J+1,1),
     $                 1, DWORK(J), DCONJG( ZWORK(J) ) )
  250    CONTINUE
C
C        Do an unshifted periodic ZQ step.
C
         DO 260  J = ILAST, JDEF+1, -1
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J,J-1,1), CS, SN, A(J,J,1) )
            A(J,J-1,1) = CZERO
            CALL ZROT( J-IFRSTM, A(IFRSTM,J,1), 1,
     $                 A(IFRSTM,J-1,1), 1, CS, SN )
            DWORK(J) = CS
            ZWORK(J) = -SN
  260    CONTINUE
         IF ( WANTQ ) THEN
            DO 270  J = ILAST, JDEF+1, -1
               CALL ZROT( N, Q(1,J-1,2), 1, Q(1,J,2),
     $                    1, DWORK(J), DCONJG( ZWORK(J) ) )
  270       CONTINUE
         END IF

C
C        Propagate the transformations through the triangular matrices.
C
         DO 310  L = 2, K
            IF ( L.GT.LDEF ) THEN
               NTRA = JDEF+2
            ELSE
               NTRA = JDEF+1
            END IF
            IF ( S(L).EQ.-1 ) THEN
               DO 280  J = ILAST, NTRA, -1
                  CS = DWORK(J)
                  SN = ZWORK(J)
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J-1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J-1,J-1,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN, A(J-1,J-1,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( ILASTM-J+1, A(J-1,J,L), LDA1, A(J,J,L),
     $                       LDA1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = SN
  280          CONTINUE
            ELSE
               DO 290  J = ILAST, NTRA, -1
                  CALL ZROT( ILASTM-J+2, A(J-1,J-1,L), LDA1,
     $                       A(J,J-1,L), LDA1, DWORK(J), ZWORK(J) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN, A(J,J,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( J-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J-1,L), 1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = -SN
  290          CONTINUE
            END IF
            IF ( WANTQ ) THEN
               LN = L+1
               IF ( L.EQ.K )  LN = 1
               DO 300  J = ILAST, NTRA, -1
                  CALL ZROT( N, Q(1,J-1,LN), 1, Q(1,J,LN), 1, DWORK(J),
     $                       DCONJG( ZWORK(J) ) )
  300          CONTINUE
            END IF
  310    CONTINUE
C
C        Apply the transformations to the left hand side of the
C        Hessenberg factor.
C
         DO 320  J = ILAST, JDEF+2, -1
            CALL ZROT( ILASTM-J+2, A(J-1,J-1,1), LDA1, A(J,J-1,1),
     $                 LDA1, DWORK(J), ZWORK(J) )
  320    CONTINUE
C
C        No QZ iteration.
C
         GOTO 450
C
C        Case III: Deflation occurs in a triangular matrix with
C                  index -1.
C
  330    CONTINUE
         JDEF = J
         IF ( JDEF.GT.( (ILAST-JLO+1)/2 ) ) THEN
C
C           Chase the zero downwards to the last position
C
            DO 350  J1 = JDEF, ILAST-1
               J = J1
               TEMP = A(J,J+1,LDEF)
               CALL ZLARTG( TEMP, A(J+1,J+1,LDEF), CS, SN,
     $                      A(J,J+1,LDEF) )
               A(J+1,J+1,LDEF) = CZERO
               CALL ZROT( ILASTM-J-1, A(J,J+2,LDEF), LDA1,
     $                    A(J+1,J+2,LDEF), LDA1, CS, SN )
               LN = LDEF+1
               IF ( LDEF.EQ.K )  LN = 1
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,LN), 1, Q(1,J+1,LN), 1, CS,
     $                       DCONJG( SN ) )
	         END IF
               DO 340  L = 1, K-1
                  IF ( LN.EQ.1 ) THEN
                     CALL ZROT( ILASTM-J+2, A(J,J-1,LN), LDA1,
     $                          A(J+1,J-1,LN), LDA1, CS, SN )
                     TEMP = A(J+1,J,LN)
                     CALL ZLARTG( TEMP, A(J+1,J-1,LN), CS, SN,
     $                            A(J+1,J,LN) )
                     A(J+1,J-1,LN) = CZERO
                     CALL ZROT( J-IFRSTM+1, A(IFRSTM,J,LN), 1,
     $                          A(IFRSTM,J-1,LN), 1, CS, SN )
	               SN = -SN
                     J = J - 1
                  ELSE IF ( S(LN).EQ.1 ) THEN
                     CALL ZROT( ILASTM-J+1, A(J,J,LN), LDA1,
     $                          A(J+1,J,LN), LDA1, CS, SN )
                     TEMP = A(J+1,J+1,LN)
                     CALL ZLARTG( TEMP, A(J+1,J,LN), CS, SN,
     $                            A(J+1,J+1,LN) )
                     A(J+1,J,LN) = CZERO
                     CALL ZROT( J-IFRSTM+1, A(IFRSTM,J+1,LN), 1,
     $                          A(IFRSTM,J,LN), 1, CS, SN )
                     SN = -SN
                  ELSE
                     CALL ZROT( J-IFRSTM+2, A(IFRSTM,J,LN), 1,
     $                          A(IFRSTM,J+1,LN), 1, CS, DCONJG( SN ) )
                     TEMP = A(J,J,LN)
                     CALL ZLARTG( TEMP, A(J+1,J,LN), CS, SN, A(J,J,LN) )
                     A(J+1,J,LN) = CZERO
                     CALL ZROT( ILASTM-J, A(J,J+1,LN), LDA1,
     $                          A(J+1,J+1,LN), LDA1, CS, SN )
                  END IF
                  LN = LN+1
                  IF ( LN.GT.K )  LN = 1
                  IF ( WANTQ ) THEN
                     CALL ZROT( N, Q(1,J,LN), 1, Q(1,J+1,LN), 1, CS,
     $                          DCONJG( SN ) )
	            END IF
  340          CONTINUE
               CALL ZROT( J-IFRSTM+1, A(IFRSTM,J,LDEF), 1,
     $                    A(IFRSTM,J+1,LDEF), 1, CS, DCONJG( SN ) )
  350       CONTINUE
C
C           Deflate the last element in the Hessenberg matrix.
C
            J = ILAST
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J,J-1,1), CS, SN, A(J,J,1) )
            A(J,J-1,1) = CZERO
            CALL ZROT( J-IFRSTM, A(IFRSTM,J,1), 1,
     $                 A(IFRSTM,J-1,1), 1, CS, SN )
            SN = -SN
   !         IF ( WANTQ.NE.0 ) THEN
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J-1,2), 1, Q(1,J,2), 1, CS,
     $                    DCONJG( SN ) )
            END IF
            DO 360  L = 2, LDEF-1
               IF ( S(L).EQ.-1 ) THEN
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J-1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J-1,J-1,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN,
     $                         A(J-1,J-1,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( ILASTM-J+1, A(J-1,J,L), LDA1,
     $                       A(J,J,L), LDA1, CS, SN )
               ELSE
                  CALL ZROT( ILASTM-J+2, A(J-1,J-1,L), LDA1,
     $                       A(J,J-1,L), LDA1, CS, SN )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN,
     $                         A(J,J,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( J-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J-1,L), 1, CS, SN )
                  SN = -SN
               END IF
               IF ( WANTQ ) THEN
                  LN = L+1
                  IF ( L.EQ.K )  LN = 1
                  CALL ZROT( N, Q(1,J-1,LN), 1, Q(1,J,LN), 1, CS,
     $                       DCONJG( SN ) )
               END IF
 360        CONTINUE
            CALL ZROT( J+1-IFRSTM, A(IFRSTM,J-1,LDEF), 1,
     $                 A(IFRSTM,J,LDEF), 1, CS, DCONJG( SN ) )
         ELSE
C
C           Chase the zero upwards to the first position.
C
            DO 380  J1 = JDEF, JLO+1,-1
               J = J1
               TEMP = A(J-1,J,LDEF)
               CALL ZLARTG( TEMP, A(J-1,J-1,LDEF), CS, SN,
     $                      A(J-1,J,LDEF) )
               A(J-1,J-1,LDEF) = CZERO
               CALL ZROT( J-IFRSTM-1, A(IFRSTM,J,LDEF), 1,
     $                    A(IFRSTM,J-1,LDEF), 1, CS, SN )
               SN = -SN
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J-1,LDEF), 1, Q(1,J,LDEF), 1, CS,
     $                       DCONJG( SN ) )
	         END IF
               LN = LDEF - 1
               DO 370  L = 1, K-1
                  IF ( LN.EQ.1 ) THEN
                     CALL ZROT( J-IFRSTM+2, A(IFRSTM,J-1,LN), 1,
     $                          A(IFRSTM,J,LN), 1, CS, DCONJG( SN ) )
                     TEMP = A(J,J-1,LN)
                     CALL ZLARTG( TEMP, A(J+1,J-1,LN), CS, SN,
     $                            A(J,J-1,LN) )
                     A(J+1,J-1,LN) = CZERO
                     CALL ZROT( ILASTM-J+1, A(J,J,LN), LDA1,
     $                          A(J+1,J,LN), LDA1, CS, SN )
                     J = J + 1
                  ELSE IF ( S(LN).EQ.-1 ) THEN
                     CALL ZROT( ILASTM-J+2, A(J-1,J-1,LN), LDA1,
     $                          A(J,J-1,LN), LDA1, CS, SN )
                     TEMP = A(J,J,LN)
                     CALL ZLARTG( TEMP, A(J,J-1,LN), CS, SN,
     $                            A(J,J,LN) )
                     A(J,J-1,LN) = CZERO
                     CALL ZROT( J-IFRSTM, A(IFRSTM,J,LN), 1,
     $                          A(IFRSTM,J-1,LN), 1, CS, SN )
                     SN = -SN
                  ELSE
                     CALL ZROT( J-IFRSTM+1, A(IFRSTM,J-1,LN), 1,
     $                          A(IFRSTM,J,LN), 1, CS, DCONJG( SN ) )
                     TEMP = A(J-1,J-1,LN)
                     CALL ZLARTG( TEMP, A(J,J-1,LN), CS, SN,
     $                            A(J-1,J-1,LN) )
                     A(J,J-1,LN) = CZERO
                     CALL ZROT( ILASTM-J+1, A(J-1,J,LN), LDA1,
     $                          A(J,J,LN), LDA1, CS, SN )
                  END IF
	            IF ( WANTQ ) THEN
                     CALL ZROT( N, Q(1,J-1,LN), 1, Q(1,J,LN), 1, CS,
     $                          DCONJG( SN ) )
	            END IF
                  LN = LN - 1
                  IF ( LN.LE.0 )  LN = K
  370          CONTINUE
               CALL ZROT( ILASTM-J+1, A(J-1,J,LDEF), LDA1, A(J,J,LDEF),
     $                    LDA1, CS, SN )
  380       CONTINUE
C
C           Deflate the first element in the Hessenberg matrix.
C
            J = JLO
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J+1,J,1), CS, SN, A(J,J,1) )
            A(J+1,J,1) = CZERO
            CALL ZROT( ILASTM-J, A(J,J+1,1), LDA1, A(J+1,J+1,1),
     $                 LDA1, CS, SN )
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, CS,
     $                    DCONJG( SN ) )
            END IF
            DO 390  L = K, LDEF+1, -1
               IF ( S(L).EQ.1 ) THEN
                  CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J+1,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                       A(J+1,J+1,L), LDA1, CS, SN )
               ELSE
                  CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1,
     $                       A(J+1,J,L), LDA1, CS, SN )
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN,
     $                         A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, SN )
                  SN = -SN
               END IF
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS,
     $                       DCONJG( SN ) )
               END IF
  390       CONTINUE
            CALL ZROT( ILASTM-J, A(J,J+1,LDEF), LDA1, A(J+1,J+1,LDEF),
     $                 LDA1, CS, SN )
         END IF
C
C        No QZ iteration.
C
         GOTO 450
C
C        Special case: A 1x1 block splits off at the bottom
C
  400    CONTINUE
         CALL ZLAPR1( BASE, K, S, A(ILAST,ILAST,1), LDA1*LDA2,
     $                ALPHA(ILAST), BETA(ILAST), SCAL(ILAST) )
C
C        Go to next block - exit if finished.
C
         ILAST = ILAST - 1
         IF ( ILAST.LT.ILO )  GOTO 470
C
C        Reset iteration counters.
C
         IITER = 0
         IF ( ZITER.NE.-1 )  ZITER = 0
         IF ( .NOT.LSCHR ) THEN
            ILASTM = ILAST
            IF ( IFRSTM.GT.ILAST )  IFRSTM = ILO
         END IF
C
C        No QZ iteration.
C
         GOTO 450        
C  
C        **************************************************************
C        *                      PERIODIC QZ STEP                      *
C        **************************************************************
C
C        It is assumed that IFIRST < ILAST.
C
  410    CONTINUE
C
         IITER = IITER + 1
         ZITER = ZITER + 1
         IF( .NOT.LSCHR ) THEN
            IFRSTM = IFIRST
         END IF
C
C        Complex single shift.
C
         IF ( ( IITER / 10 )*10.EQ.IITER ) THEN
C
C           Exceptional shift.
C
            CALL ZLARNV( 2, ISEED, 2, RND )
            CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
         ELSE
	      CALL ZLARTG( CONE, CONE, CS, SN, TEMP )
            DO 420  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
	            CALL ZLARTG( A(IFIRST,IFIRST,L)*CS,
     $                         A(ILAST,ILAST,L)*DCONJG(SN),
     $                         CS, SN, TEMP )
               ELSE
	            CALL ZLARTG( A(ILAST,ILAST,L)*CS,
     $                         -A(IFIRST,IFIRST,L)*DCONJG(SN),
     $                         CS, SN, TEMP )
	            SN = -SN
               END IF
  420       CONTINUE
            CALL ZLARTG( A(IFIRST,IFIRST,1)*CS
     $                   -DCONJG(SN)*A(ILAST,ILAST,1),
     $                   A(IFIRST+1,IFIRST,1)*CS, CS, SN, TEMP )
         END IF
C
C        Do the sweeps.
C
         DO 440  J1 = IFIRST-1, ILAST-2
            J = J1 + 1
C
C           Create bulge if J1 = IFIRST - 1, otherwise chase bulge.
C
            IF ( J1.LT.IFIRST ) THEN
               CALL ZROT( ILASTM-J+1, A(J,J,1), LDA1, A(J+1,J,1), LDA1,
     $                    CS, SN )
            ELSE
               TEMP = A(J,J-1,1)
               CALL ZLARTG( TEMP, A(J+1,J-1,1), CS, SN, A(J,J-1,1) )
               A(J+1,J-1,1) = CZERO
               CALL ZROT( ILASTM-J+1, A(J,J,1), LDA1, A(J+1,J,1),
     $                    LDA1, CS, SN )
            END IF
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, CS,
     $                    DCONJG( SN ) )
            END IF
C
C           Propagate rotation through AK, ..., A2 to A1.
C
            DO 430  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
                  CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J+1,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                       A(J+1,J+1,L), LDA1, CS, SN )
               ELSE
                  CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1, A(J+1,J,L),
     $                       LDA1, CS, SN )
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, SN )
                  SN = -SN
               END IF
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS,
     $                       DCONJG( SN ) )
               END IF
  430       CONTINUE
            CALL ZROT( MIN(J+2,ILASTM)-IFRSTM+1, A(IFRSTM,J,1), 1,
     $                 A(IFRSTM,J+1,1), 1, CS, DCONJG( SN ) )
  440    CONTINUE
C
C        End of iteration loop.
C
  450    CONTINUE
  460 CONTINUE
C
C     Drop through = non-convergence
C
      INFO = ILAST
      GO TO 550
C
C     Successful completion of all QZ steps
C
  470 CONTINUE
C
C     Set eigenvalues 1:ILO-1
C
      DO 480  J = 1,ILO-1
         CALL ZLAPR1( BASE, K, S, A(J,J,1), LDA1*LDA2, ALPHA(J),
     $                BETA(J), SCAL(J) )
  480 CONTINUE
      IF ( LSCHR ) THEN
C
C        Scale A(2,:,:) .. A(K,:,:).
C
         DO 540  L = K, 2, -1
	      IF ( S(L).EQ.1 )  THEN
	         DO 490 J = 1, N
                  ABST = ABS( A(J,J,L) )
                  IF ( ABST.GT.SAFMIN ) THEN
                     TEMP = DCONJG( A(J,J,L) / ABST )
                     A(J,J,L ) = ABST
                     IF ( J.LT.N )
     $                  CALL ZSCAL( N-J, TEMP, A(J,J+1,L), LDA1 )
	            ELSE
	               TEMP = CONE
	            END IF
                  ZWORK(J) = TEMP
  490          CONTINUE	
	      ELSE
               DO 500  J = 1, N
                  ABST = ABS( A(J,J,L) )
                  IF ( ABST.GT.SAFMIN ) THEN
                     TEMP = DCONJG( A(J,J,L) / ABST )
                     A(J,J,L ) = ABST
                     CALL ZSCAL( J-1, TEMP, A(1,J,L), 1 )
	            ELSE
	               TEMP = CONE
	            END IF
                  ZWORK(J) = DCONJG(TEMP)
  500          CONTINUE	
	      END IF
            IF ( WANTQ ) THEN
	         DO 510  J = 1, N
                  CALL ZSCAL( N, DCONJG( ZWORK(J) ), Q(1,J,L), 1 )
  510          CONTINUE
            END IF
         	IF ( S(L-1).EQ.1 )  THEN
	         DO 520  J = 1, N
	            CALL ZSCAL( J, DCONJG( ZWORK(J) ), A(1,J,L-1), 1 )
  520          CONTINUE
            ELSE
	         DO 530  J = 1, N
	            CALL ZSCAL( N-J+1, ZWORK(J), A(J,J,L-1), LDA1 )
  530          CONTINUE              
	      END IF
  540    CONTINUE
	END IF
      INFO = 0
C
  550 CONTINUE
C
      DWORK(1) = DBLE( N )
      ZWORK(1) = DCMPLX( N, 0 )
      RETURN     
C *** Last line of ZPGEQZ ***
      END
      SUBROUTINE ZPGORD( WANTQ, K, N, S, SELECT, A, LDA1, LDA2,
     $                   ALPHA, BETA, SCAL, Q, LDQ1, LDQ2, M,
     $                   ZWORK, LZWORK, INFO )
	IMPLICIT NONE
C
C     PURPOSE
C
C     ZPGORD reorders the periodic Schur decomposition of a complex
C     generalized matrix product
C
C                                 S(2)                 S(K)
C          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K),
C
C     (in terms of unitary equivalence transformations), so that a
C     selected cluster of eigenvalues appears in the leading diagonal
C     blocks of the matrix product. The leading columns of the
C     orthogonal factors contained in Q form unitary bases of the
C     corresponding periodic eigenspaces (deflating subspaces).
C     A must be in periodic Schur form, that is, all factors of A must
C     be upper triangular.
C
C     If WANTQ = .TRUE., then the unitary factors are computed and
C     stored in the array Q so that for S(I) = 1,
C
C                         H
C             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)
C                          H                                        (1)
C         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),
C
C     and for S(I) = -1,
C
C                                  H
C             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)
C                                   H                               (2)
C         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).
C
C     ARGUMEMTS
C
C     Mode Parameters
C
C     WANTQ   (input) LOGICAL
C             = .FALSE.: do not modify Q;
C             = .TRUE. : modify the array Q by the unitary
C                        transformations that are applied to the
C                        matrices in A for reordering.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of factors.  K >= 1.
C
C     N       (input)  INTEGER
C             The order of each factor in A.  N >= 0.
C
C     S       (input)  INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     SELECT  (input) LOGICAL array, dimension (N)
C             SELECT specifies the eigenvalues in the selected cluster.
C             To select the eigenvalue corresponding to the (j,j)
C             diagonal entries, SELECT(j) must be set to .TRUE..
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array
C             must contain the factors in periodic Schur form
C             form, that is, all factors are upper triangular.
C             On exit, if INFO = 0, the leading N-by-N-by-K part of
C             this array contains the factors of the reordered periodic
C             Schur form. Moreover, A(:,:,2), ..., A(:,:,K) are
C             normalized so that their diagonals contain nonnegative
C             real numbers.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A.  LDA1 >= MAX(1,N).
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A.  LDA2 >= MAX(1,N).
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             On entry, if WANTQ = .TRUE., the leading N-by-N-by-K part
C             of this array must contain the initial unitary factors
C             as described in (1)-(2).
C             On exit, if WANTQ = .TRUE., the leading N-by-N-by-K part
C             of this array contains the modified orthogonal factors as
C             described in (1)-(2).
C
C     LDQ1    (input)  INTEGER
C             The first leading dimension of Q.
C             If WANTQ = .TRUE.,  LDQ1 >= MAX(1,N).
C
C     LDQ2    (input)  INTEGER
C             The second leading dimension of Q.
C             If WANTQ = .TRUE.,  LDQ2 >= MAX(1,N).
C
C     M       (output) INTEGER
C             The dimension of the specified periodic eigenspace.
C
C     Workspace
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the minimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.  LZWORK >= MAX( K, N ).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0       : succesful exit;
C             < 0       : if INFO = -i, the i-th argument had an
C                         illegal value;
C             = 1       : the periodic QZ algorithm failed to converge.
C
C     METHOD
C
C     A complex version of the periodic QZ algorithm [1] with 
C     perfect shifts is used. For more details see [2]. This is not
C     a safe method. It is advisable to check whether the eigenvalues
C     are really reorderd.
C
C     REFERENCES
C
C     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P.
C         The periodic Schur decomposition; algorithm and applications.
C         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992.
C
C     [2] Kressner, D.
C         An efficient and reliable implementation of the periodic QZ
C         algorithm. In IFAC Workshop on Periodic Control Systems, 2001.
C
C     NUMERICAL ASPECTS
C
C     The implemented method is numerically backward stable.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      LOGICAL           WANTQ
      INTEGER           K, INFO, LDA1, LDA2, LDQ1, LDQ2, LZWORK, M, N
C     .. Array Arguments ..
      LOGICAL           SELECT(*)
      INTEGER           S(*), SCAL(*)
      COMPLEX*16        A(LDA1, LDA2, *), ALPHA(*), BETA(*),
     $                  Q(LDQ1, LDQ2, *), ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           SOK
      INTEGER           I, IERR, J, JS, L
	DOUBLE PRECISION  ABST, BASE, SAFMIN
	COMPLEX*16        TEMP
C     .. Local Arrays ..
	DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          ZLAPR1, ZPGEX2, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, DCONJG, MAX
C
C     .. Executable Statements ..
C
      INFO = 0
	SAFMIN = DLAMCH( 'SafeMinimum' )
      BASE = DLAMCH( 'Base' )
C
C     Check the scalar input parameters.
C
      IF ( K.LT.1 ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE
	   SOK = S(1).EQ.1
	   DO 10  L = 2, K
	      SOK = S(L).EQ.1 .OR. S(L).EQ.-1
   10    CONTINUE
         IF ( .NOT.SOK ) THEN
	      INFO = -4
         ELSE IF ( LDA1 .LT. MAX( 1, N ) ) THEN
            INFO = -7
         ELSE IF ( LDA2 .LT. MAX( 1, N ) ) THEN	
            INFO = -8
         ELSE IF ( WANTQ .AND. LDQ1.LT.MAX( 1, N ) ) THEN
            INFO = -13
         ELSE IF ( WANTQ .AND. LDQ2.LT.MAX( 1, N ) ) THEN
            INFO = -14
	   ELSE IF ( LZWORK.LT.MAX( K, N ) ) THEN
            INFO = -17
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPGORD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         ZWORK(1) = CONE
         RETURN
      END IF
C
C     Set M to the dimension of the specified deflating subspace.
C
      M = 0
      DO 20 J = 1, N
         IF ( SELECT( J ) )
     $      M = M + 1
   20 CONTINUE
C
      JS = 0
      DO 40 J = 1, N
         IF ( SELECT( J ) ) THEN
            JS = JS + 1
C
C           Swap the J-th block to position JS.
C
            IF ( J.NE.JS ) THEN
	         DO 30  I = J-1, JS, -1
                  CALL ZPGEX2( .TRUE., K, N, I, S, A, LDA1, LDA2, Q,
     $                         LDQ1, LDQ2, ZWORK, IERR )
	            IF ( IERR.NE.0 ) THEN
	               INFO = 1
	               GO TO 120
	            END IF
   30          CONTINUE
	      END IF
         END IF
   40 CONTINUE
C
C     Rescale matrices and recompute eigenvalues.
C
      DO 100  L = K, 2, -1
	   IF ( S(L).EQ.1 )  THEN
	      DO 50 J = 1, N
               ABST = ABS( A(J,J,L) )
               IF ( ABST.GT.SAFMIN ) THEN
                  TEMP = DCONJG( A(J,J,L) / ABST )
                  A(J,J,L ) = ABST
                  CALL ZSCAL( N-J, TEMP, A(J,J+1,L), LDA1 )
	         ELSE
	            TEMP = CONE
	         END IF
               ZWORK(J) = TEMP
   50       CONTINUE	
	   ELSE
            DO 60  J = 1, N
               ABST = ABS( A(J,J,L) )
               IF ( ABST.GT.SAFMIN ) THEN
                  TEMP = DCONJG( A(J,J,L) / ABST )
                  A(J,J,L ) = ABST
                  CALL ZSCAL( J-1, TEMP, A(1,J,L), 1 )
	         ELSE
	            TEMP = CONE
	         END IF
               ZWORK(J) = DCONJG(TEMP)
   60       CONTINUE	
	   END IF
         IF ( WANTQ ) THEN
	      DO 70  J = 1, N
               CALL ZSCAL( N, DCONJG( ZWORK(J) ), Q(1,J,L), 1 )
   70       CONTINUE
         END IF
         IF ( S(L-1).EQ.1 )  THEN
	      DO 80  J = 1, N
	         CALL ZSCAL( J, DCONJG( ZWORK(J) ), A(1,J,L-1), 1 )
   80       CONTINUE
         ELSE
	      DO 90  J = 1, N
	         CALL ZSCAL( N-J+1, ZWORK(J), A(J,J,L-1), LDA1 )
   90       CONTINUE              
	   END IF
  100 CONTINUE
      DO 110  J = 1, N
         CALL ZLAPR1( BASE, K, S, A(J,J,1), LDA1*LDA2, ALPHA(J),
     $                BETA(J), SCAL(J) )
  110 CONTINUE
C
  120 CONTINUE
      ZWORK(1) = DCMPLX( MAX( K, N ), 0 )
      RETURN
C *** Last line of ZPGORD ***
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     Modified to be used from Matlab mexfiles (February, 1999).
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      CHARACTER*72 ERRTXT
      WRITE( ERRTXT, FMT = 9999 ) SRNAME, INFO
*
      CALL mexErrMsgTxt( ERRTXT )
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
      RETURN
*
*     End of XERBLA
*
      END
