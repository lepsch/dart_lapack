      SUBROUTINE SGTT01( N, DL, D, DU, DLF, DF, DUF, DU2, IPIV, WORK,
     $                   LDWORK, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDWORK, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               D( * ), DF( * ), DL( * ), DLF( * ), DU( * ),
     $                   DU2( * ), DUF( * ), RWORK( * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IP, J, LASTJ
      REAL               ANORM, EPS, LI
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGT, SLANHS
      EXTERNAL           SLAMCH, SLANGT, SLANHS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SSWAP
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      EPS = SLAMCH( 'Epsilon' )
*
*     Copy the matrix U to WORK.
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            WORK( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, N
         IF( I.EQ.1 ) THEN
            WORK( I, I ) = DF( I )
            IF( N.GE.2 )
     $         WORK( I, I+1 ) = DUF( I )
            IF( N.GE.3 )
     $         WORK( I, I+2 ) = DU2( I )
         ELSE IF( I.EQ.N ) THEN
            WORK( I, I ) = DF( I )
         ELSE
            WORK( I, I ) = DF( I )
            WORK( I, I+1 ) = DUF( I )
            IF( I.LT.N-1 )
     $         WORK( I, I+2 ) = DU2( I )
         END IF
   30 CONTINUE
*
*     Multiply on the left by L.
*
      LASTJ = N
      DO 40 I = N - 1, 1, -1
         LI = DLF( I )
         CALL SAXPY( LASTJ-I+1, LI, WORK( I, I ), LDWORK,
     $               WORK( I+1, I ), LDWORK )
         IP = IPIV( I )
         IF( IP.EQ.I ) THEN
            LASTJ = MIN( I+2, N )
         ELSE
            CALL SSWAP( LASTJ-I+1, WORK( I, I ), LDWORK, WORK( I+1, I ),
     $                  LDWORK )
         END IF
   40 CONTINUE
*
*     Subtract the matrix A.
*
      WORK( 1, 1 ) = WORK( 1, 1 ) - D( 1 )
      IF( N.GT.1 ) THEN
         WORK( 1, 2 ) = WORK( 1, 2 ) - DU( 1 )
         WORK( N, N-1 ) = WORK( N, N-1 ) - DL( N-1 )
         WORK( N, N ) = WORK( N, N ) - D( N )
         DO 50 I = 2, N - 1
            WORK( I, I-1 ) = WORK( I, I-1 ) - DL( I-1 )
            WORK( I, I ) = WORK( I, I ) - D( I )
            WORK( I, I+1 ) = WORK( I, I+1 ) - DU( I )
   50    CONTINUE
      END IF
*
*     Compute the 1-norm of the tridiagonal matrix A.
*
      ANORM = SLANGT( '1', N, DL, D, DU )
*
*     Compute the 1-norm of WORK, which is only guaranteed to be
*     upper Hessenberg.
*
      RESID = SLANHS( '1', N, WORK, LDWORK, RWORK )
*
*     Compute norm(L*U - A) / (norm(A) * EPS)
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         RESID = ( RESID / ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of SGTT01
*
      END