      SUBROUTINE CLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IJOB, LDZ, N;
      REAL               RDSCAL, RDSUM
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      COMPLEX            RHS( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                MAXDIM;
      PARAMETER          ( MAXDIM = 2 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      int                I, INFO, J, K;
      REAL               RTEMP, SCALE, SMINU, SPLUS
      COMPLEX            BM, BP, PMONE, TEMP
*     ..
*     .. Local Arrays ..
      REAL               RWORK( MAXDIM )
      COMPLEX            WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CGECON, CGESC2, CLASSQ, CLASWP, CSCAL
*     ..
*     .. External Functions ..
      REAL               SCASUM
      COMPLEX            CDOTC
      EXTERNAL           SCASUM, CDOTC
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( IJOB.NE.2 ) THEN
*
*        Apply permutations IPIV to RHS
*
         CALL CLASWP( 1, RHS, LDZ, 1, N-1, IPIV, 1 )
*
*        Solve for L-part choosing RHS either to +1 or -1.
*
         PMONE = -CONE
         DO 10 J = 1, N - 1
            BP = RHS( J ) + CONE
            BM = RHS( J ) - CONE
            SPLUS = ONE
*
*           Look-ahead for L- part RHS(1:N-1) = +-1
*           SPLUS and SMIN computed more efficiently than in BSOLVE[1].
*
            SPLUS = SPLUS + REAL( CDOTC( N-J, Z( J+1, J ), 1, Z( J+1, J ), 1 ) )
            SMINU = REAL( CDOTC( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 ) )
            SPLUS = SPLUS*REAL( RHS( J ) )
            IF( SPLUS.GT.SMINU ) THEN
               RHS( J ) = BP
            ELSE IF( SMINU.GT.SPLUS ) THEN
               RHS( J ) = BM
            ELSE
*
*              In this case the updating sums are equal and we can
*              choose RHS(J) +1 or -1. The first time this happens we
*              choose -1, thereafter +1. This is a simple way to get
*              good estimates of matrices like Byers well-known example
*              (see [1]). (Not done in BSOLVE.)
*
               RHS( J ) = RHS( J ) + PMONE
               PMONE = CONE
            END IF
*
*           Compute the remaining r.h.s.
*
            TEMP = -RHS( J )
            CALL CAXPY( N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 )
   10    CONTINUE
*
*        Solve for U- part, lockahead for RHS(N) = +-1. This is not done
*        In BSOLVE and will hopefully give us a better estimate because
*        any ill-conditioning of the original matrix is transferred to U
*        and not to L. U(N, N) is an approximation to sigma_min(LU).
*
         CALL CCOPY( N-1, RHS, 1, WORK, 1 )
         WORK( N ) = RHS( N ) + CONE
         RHS( N ) = RHS( N ) - CONE
         SPLUS = ZERO
         SMINU = ZERO
         DO 30 I = N, 1, -1
            TEMP = CONE / Z( I, I )
            WORK( I ) = WORK( I )*TEMP
            RHS( I ) = RHS( I )*TEMP
            DO 20 K = I + 1, N
               WORK( I ) = WORK( I ) - WORK( K )*( Z( I, K )*TEMP )
               RHS( I ) = RHS( I ) - RHS( K )*( Z( I, K )*TEMP )
   20       CONTINUE
            SPLUS = SPLUS + ABS( WORK( I ) )
            SMINU = SMINU + ABS( RHS( I ) )
   30    CONTINUE
         IF( SPLUS.GT.SMINU ) CALL CCOPY( N, WORK, 1, RHS, 1 )
*
*        Apply the permutations JPIV to the computed solution (RHS)
*
         CALL CLASWP( 1, RHS, LDZ, 1, N-1, JPIV, -1 )
*
*        Compute the sum of squares
*
         CALL CLASSQ( N, RHS, 1, RDSCAL, RDSUM )
         RETURN
      END IF
*
*     ENTRY IJOB = 2
*
*     Compute approximate nullvector XM of Z
*
      CALL CGECON( 'I', N, Z, LDZ, ONE, RTEMP, WORK, RWORK, INFO )
      CALL CCOPY( N, WORK( N+1 ), 1, XM, 1 )
*
*     Compute RHS
*
      CALL CLASWP( 1, XM, LDZ, 1, N-1, IPIV, -1 )
      TEMP = CONE / SQRT( CDOTC( N, XM, 1, XM, 1 ) )
      CALL CSCAL( N, TEMP, XM, 1 )
      CALL CCOPY( N, XM, 1, XP, 1 )
      CALL CAXPY( N, CONE, RHS, 1, XP, 1 )
      CALL CAXPY( N, -CONE, XM, 1, RHS, 1 )
      CALL CGESC2( N, Z, LDZ, RHS, IPIV, JPIV, SCALE )
      CALL CGESC2( N, Z, LDZ, XP, IPIV, JPIV, SCALE )
      IF( SCASUM( N, XP, 1 ).GT.SCASUM( N, RHS, 1 ) ) CALL CCOPY( N, XP, 1, RHS, 1 )
*
*     Compute the sum of squares
*
      CALL CLASSQ( N, RHS, 1, RDSCAL, RDSUM )
      RETURN
*
*     End of CLATDF
*
      END
