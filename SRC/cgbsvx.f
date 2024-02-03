      SUBROUTINE CGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*  Moved setting of INFO = N+1 so INFO does not subsequently get
*  overwritten.  Sven, 17 Mar 05.
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU
      CHARACTER          NORM
      INTEGER            I, INFEQU, J, J1, J2
      REAL               AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, RPVGRW, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGB, CLANTB, SLAMCH
      EXTERNAL           LSAME, CLANGB, CLANTB, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGBCON, CGBEQU, CGBRFS, CGBTRF, CGBTRS, CLACPY, CLAQGB, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      ELSE
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF
*
*     Test the input parameters.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KL.LT.0 ) THEN
         INFO = -4
      ELSE IF( KU.LT.0 ) THEN
         INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -8
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
         INFO = -10
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT. ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -12
      ELSE
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
   10       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -13
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
   20       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -14
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               COLCND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -16
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -18
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGBSVX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL CGBEQU( N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL CLAQGB( N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
      END IF
*
*     Scale the right hand side.
*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) THEN
            DO 40 J = 1, NRHS
               DO 30 I = 1, N
                  B( I, J ) = R( I )*B( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( COLEQU ) THEN
         DO 60 J = 1, NRHS
            DO 50 I = 1, N
               B( I, J ) = C( I )*B( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the LU factorization of the band matrix A.
*
         DO 70 J = 1, N
            J1 = MAX( J-KU, 1 )
            J2 = MIN( J+KL, N )
            CALL CCOPY( J2-J1+1, AB( KU+1-J+J1, J ), 1, AFB( KL+KU+1-J+J1, J ), 1 )
   70    CONTINUE
*
         CALL CGBTRF( N, N, KL, KU, AFB, LDAFB, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.GT.0 ) THEN
*
*           Compute the reciprocal pivot growth factor of the
*           leading rank-deficient INFO columns of A.
*
            ANORM = ZERO
            DO 90 J = 1, INFO
               DO 80 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )
                  ANORM = MAX( ANORM, ABS( AB( I, J ) ) )
   80          CONTINUE
   90       CONTINUE
            RPVGRW = CLANTB( 'M', 'U', 'N', INFO, MIN( INFO-1, KL+KU ), AFB( MAX( 1, KL+KU+2-INFO ), 1 ), LDAFB, RWORK )
            IF( RPVGRW.EQ.ZERO ) THEN
               RPVGRW = ONE
            ELSE
               RPVGRW = ANORM / RPVGRW
            END IF
            RWORK( 1 ) = RPVGRW
            RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A and the
*     reciprocal pivot growth factor RPVGRW.
*
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = CLANGB( NORM, N, KL, KU, AB, LDAB, RWORK )
      RPVGRW = CLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, RWORK )
      IF( RPVGRW.EQ.ZERO ) THEN
         RPVGRW = ONE
      ELSE
         RPVGRW = CLANGB( 'M', N, KL, KU, AB, LDAB, RWORK ) / RPVGRW
      END IF
*
*     Compute the reciprocal of the condition number of A.
*
      CALL CGBCON( NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND, WORK, RWORK, INFO )
*
*     Compute the solution matrix X.
*
      CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL CGBTRS( TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL CGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      IF( NOTRAN ) THEN
         IF( COLEQU ) THEN
            DO 110 J = 1, NRHS
               DO 100 I = 1, N
                  X( I, J ) = C( I )*X( I, J )
  100          CONTINUE
  110       CONTINUE
            DO 120 J = 1, NRHS
               FERR( J ) = FERR( J ) / COLCND
  120       CONTINUE
         END IF
      ELSE IF( ROWEQU ) THEN
         DO 140 J = 1, NRHS
            DO 130 I = 1, N
               X( I, J ) = R( I )*X( I, J )
  130       CONTINUE
  140    CONTINUE
         DO 150 J = 1, NRHS
            FERR( J ) = FERR( J ) / ROWCND
  150    CONTINUE
      END IF
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RWORK( 1 ) = RPVGRW
      RETURN
*
*     End of CGBSVX
*
      END
