      SUBROUTINE DLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, ISEED, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS, UPLO, XTYPE
      CHARACTER*3        PATH
      int                INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI
      CHARACTER          C1, DIAG
      CHARACTER*2        C2
      int                J, MB, NX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           LSAME, LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGBMV, DGEMM, DLACPY, DLARNV, DSBMV, DSPMV, DSYMM, DTBMV, DTPMV, DTRMM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )
      TRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      NOTRAN = .NOT.TRAN
      GEN = LSAME( PATH( 2: 2 ), 'G' )
      QRS = LSAME( PATH( 2: 2 ), 'Q' ) .OR. LSAME( PATH( 3: 3 ), 'Q' )
      SYM = LSAME( PATH( 2: 2 ), 'P' ) .OR. LSAME( PATH( 2: 2 ), 'S' )
      TRI = LSAME( PATH( 2: 2 ), 'T' )
      BAND = LSAME( PATH( 3: 3 ), 'B' )
      IF( .NOT.LSAME( C1, 'Double precision' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( XTYPE, 'N' ) .OR. LSAME( XTYPE, 'C' ) ) ) THEN
         INFO = -2
      ELSE IF( ( SYM .OR. TRI ) .AND. .NOT. ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( ( GEN .OR. QRS ) .AND. .NOT. ( TRAN .OR. LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( BAND .AND. KL.LT.0 ) THEN
         INFO = -7
      ELSE IF( BAND .AND. KU.LT.0 ) THEN
         INFO = -8
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -9
      ELSE IF( ( .NOT.BAND .AND. LDA.LT.MAX( 1, M ) ) .OR. ( BAND .AND. ( SYM .OR. TRI ) .AND. LDA.LT.KL+1 ) .OR. ( BAND .AND. GEN .AND. LDA.LT.KL+KU+1 ) ) THEN
         INFO = -11
      ELSE IF( ( NOTRAN .AND. LDX.LT.MAX( 1, N ) ) .OR. ( TRAN .AND. LDX.LT.MAX( 1, M ) ) ) THEN
         INFO = -13
      ELSE IF( ( NOTRAN .AND. LDB.LT.MAX( 1, M ) ) .OR. ( TRAN .AND. LDB.LT.MAX( 1, N ) ) ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLARHS', -INFO )
         RETURN
      END IF
*
*     Initialize X to NRHS random vectors unless XTYPE = 'C'.
*
      IF( TRAN ) THEN
         NX = M
         MB = N
      ELSE
         NX = N
         MB = M
      END IF
      IF( .NOT.LSAME( XTYPE, 'C' ) ) THEN
         DO 10 J = 1, NRHS
            CALL DLARNV( 2, ISEED, N, X( 1, J ) )
   10    CONTINUE
      END IF
*
*     Multiply X by op(A) using an appropriate
*     matrix multiply routine.
*
      IF( LSAMEN( 2, C2, 'GE' ) .OR. LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C2, 'RQ' ) ) THEN
*
*        General matrix
*
         CALL DGEMM( TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'SY' ) ) THEN
*
*        Symmetric matrix, 2-D storage
*
         CALL DSYMM( 'Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
*        General matrix, band storage
*
         DO 20 J = 1, NRHS
            CALL DGBMV( TRANS, MB, NX, KL, KU, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   20    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
*
*        Symmetric matrix, band storage
*
         DO 30 J = 1, NRHS
            CALL DSBMV( UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   30    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'PP' ) .OR. LSAMEN( 2, C2, 'SP' ) ) THEN
*
*        Symmetric matrix, packed storage
*
         DO 40 J = 1, NRHS
            CALL DSPMV( UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   40    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN
*
*        Triangular matrix.  Note that for triangular matrices,
*           KU = 1 => non-unit triangular
*           KU = 2 => unit triangular
*
         CALL DLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         IF( KU.EQ.2 ) THEN
            DIAG = 'U'
         ELSE
            DIAG = 'N'
         END IF
         CALL DTRMM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'TP' ) ) THEN
*
*        Triangular matrix, packed storage
*
         CALL DLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         IF( KU.EQ.2 ) THEN
            DIAG = 'U'
         ELSE
            DIAG = 'N'
         END IF
         DO 50 J = 1, NRHS
            CALL DTPMV( UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 )
   50    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
*
*        Triangular matrix, banded storage
*
         CALL DLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         IF( KU.EQ.2 ) THEN
            DIAG = 'U'
         ELSE
            DIAG = 'N'
         END IF
         DO 60 J = 1, NRHS
            CALL DTBMV( UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 )
   60    CONTINUE
*
      ELSE
*
*        If PATH is none of the above, return with an error code.
*
         INFO = -1
         CALL XERBLA( 'DLARHS', -INFO )
      END IF
*
      RETURN
*
*     End of DLARHS
*
      END
