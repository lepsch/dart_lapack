      double           FUNCTION DQPT01( M, N, K, A, AF, LDA, TAU, JPVT, WORK, LWORK );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      int                JPVT( * )
      double             A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      int                I, INFO, J
      double             NORMA;
*     ..
*     .. Local Arrays ..
      double             RWORK( 1 );
*     ..
*     .. External Functions ..
      double             DLAMCH, DLANGE;
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DORMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      DQPT01 = ZERO
*
*     Test if there is enough workspace
*
      IF( LWORK.LT.M*N+N ) THEN
         CALL XERBLA( 'DQPT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
      DO J = 1, K
*
*        Copy the upper triangular part of the factor R stored
*        in AF(1:K,1:K) into the work array WORK.
*
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = AF( I, J )
         END DO
*
*        Zero out the elements below the diagonal in the work array.
*
         DO I = J + 1, M
            WORK( ( J-1 )*M+I ) = ZERO
         END DO
      END DO
*
*     Copy columns (K+1,N) from AF into the work array WORK.
*     AF(1:K,K+1:N) contains the rectangular block of the upper trapezoidal
*     factor R, AF(K+1:M,K+1:N) contains the partially updated residual
*     matrix of R.
*
      DO J = K + 1, N
         CALL DCOPY( M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 )
      END DO
*
      CALL DORMQR( 'Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
      DO J = 1, N
*
*        Compare J-th column of QR and JPVT(J)-th column of A.
*
         CALL DAXPY( M, -ONE, A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 )
      END DO
*
      DQPT01 = DLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( DBLE( MAX( M, N ) )*DLAMCH( 'Epsilon' ) )       IF( NORMA.NE.ZERO ) DQPT01 = DQPT01 / NORMA
*
      RETURN
*
*     End of DQPT01
*
      END
