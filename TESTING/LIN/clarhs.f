      SUBROUTINE CLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS,
     $                   A, LDA, X, LDX, B, LDB, ISEED, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI
      CHARACTER          C1, DIAG
      CHARACTER*2        C2
      INTEGER            J, MB, NX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           LSAME, LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGBMV, CGEMM, CHBMV, CHEMM, CHPMV, CLACPY,
     $                   CLARNV, CSBMV, CSPMV, CSYMM, CTBMV, CTPMV,
     $                   CTRMM, XERBLA
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
      SYM = LSAME( PATH( 2: 2 ), 'P' ) .OR.
     $      LSAME( PATH( 2: 2 ), 'S' ) .OR. LSAME( PATH( 2: 2 ), 'H' )
      TRI = LSAME( PATH( 2: 2 ), 'T' )
      BAND = LSAME( PATH( 3: 3 ), 'B' )
      IF( .NOT.LSAME( C1, 'Complex precision' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( XTYPE, 'N' ) .OR. LSAME( XTYPE, 'C' ) ) )
     $          THEN
         INFO = -2
      ELSE IF( ( SYM .OR. TRI ) .AND. .NOT.
     $         ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( ( GEN.OR.QRS ) .AND.
     $   .NOT.( TRAN .OR. LSAME( TRANS, 'N' ) ) ) THEN
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
      ELSE IF( ( .NOT.BAND .AND. LDA.LT.MAX( 1, M ) ) .OR.
     $         ( BAND .AND. ( SYM .OR. TRI ) .AND. LDA.LT.KL+1 ) .OR.
     $         ( BAND .AND. GEN .AND. LDA.LT.KL+KU+1 ) ) THEN
         INFO = -11
      ELSE IF( ( NOTRAN .AND. LDX.LT.MAX( 1, N ) ) .OR.
     $         ( TRAN .AND. LDX.LT.MAX( 1, M ) ) ) THEN
         INFO = -13
      ELSE IF( ( NOTRAN .AND. LDB.LT.MAX( 1, M ) ) .OR.
     $         ( TRAN .AND. LDB.LT.MAX( 1, N ) ) ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLARHS', -INFO )
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
            CALL CLARNV( 2, ISEED, N, X( 1, J ) )
   10    CONTINUE
      END IF
*
*     Multiply X by op(A) using an appropriate
*     matrix multiply routine.
*
      IF( LSAMEN( 2, C2, 'GE' ) .OR. LSAMEN( 2, C2, 'QR' ) .OR.
     $    LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR.
     $    LSAMEN( 2, C2, 'RQ' ) ) THEN
*
*        General matrix
*
         CALL CGEMM( TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX,
     $               ZERO, B, LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'HE' ) ) THEN
*
*        Hermitian matrix, 2-D storage
*
         CALL CHEMM( 'Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO,
     $               B, LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN
*
*        Symmetric matrix, 2-D storage
*
         CALL CSYMM( 'Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO,
     $               B, LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
*        General matrix, band storage
*
         DO 20 J = 1, NRHS
            CALL CGBMV( TRANS, M, N, KL, KU, ONE, A, LDA, X( 1, J ), 1,
     $                  ZERO, B( 1, J ), 1 )
   20    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'PB' ) .OR. LSAMEN( 2, C2, 'HB' ) ) THEN
*
*        Hermitian matrix, band storage
*
         DO 30 J = 1, NRHS
            CALL CHBMV( UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO,
     $                  B( 1, J ), 1 )
   30    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'SB' ) ) THEN
*
*        Symmetric matrix, band storage
*
         DO 40 J = 1, NRHS
            CALL CSBMV( UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO,
     $                  B( 1, J ), 1 )
   40    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'PP' ) .OR. LSAMEN( 2, C2, 'HP' ) ) THEN
*
*        Hermitian matrix, packed storage
*
         DO 50 J = 1, NRHS
            CALL CHPMV( UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ),
     $                  1 )
   50    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'SP' ) ) THEN
*
*        Symmetric matrix, packed storage
*
         DO 60 J = 1, NRHS
            CALL CSPMV( UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ),
     $                  1 )
   60    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN
*
*        Triangular matrix.  Note that for triangular matrices,
*           KU = 1 => non-unit triangular
*           KU = 2 => unit triangular
*
         CALL CLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         IF( KU.EQ.2 ) THEN
            DIAG = 'U'
         ELSE
            DIAG = 'N'
         END IF
         CALL CTRMM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,
     $               LDB )
*
      ELSE IF( LSAMEN( 2, C2, 'TP' ) ) THEN
*
*        Triangular matrix, packed storage
*
         CALL CLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         IF( KU.EQ.2 ) THEN
            DIAG = 'U'
         ELSE
            DIAG = 'N'
         END IF
         DO 70 J = 1, NRHS
            CALL CTPMV( UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 )
   70    CONTINUE
*
      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
*
*        Triangular matrix, banded storage
*
         CALL CLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         IF( KU.EQ.2 ) THEN
            DIAG = 'U'
         ELSE
            DIAG = 'N'
         END IF
         DO 80 J = 1, NRHS
            CALL CTBMV( UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 )
   80    CONTINUE
*
      ELSE
*
*        If none of the above, set INFO = -1 and return
*
         INFO = -1
         CALL XERBLA( 'CLARHS', -INFO )
      END IF
*
      RETURN
*
*     End of CLARHS
*
      END