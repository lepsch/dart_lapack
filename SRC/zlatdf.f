      SUBROUTINE ZLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IJOB, LDZ, N;
      double             RDSCAL, RDSUM;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      COMPLEX*16         RHS( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                MAXDIM;
      const              MAXDIM = 2 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, K;
      double             RTEMP, SCALE, SMINU, SPLUS;
      COMPLEX*16         BM, BP, PMONE, TEMP
      // ..
      // .. Local Arrays ..
      double             RWORK( MAXDIM );
      COMPLEX*16         WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZGECON, ZGESC2, ZLASSQ, ZLASWP, ZSCAL
      // ..
      // .. External Functions ..
      double             DZASUM;
      COMPLEX*16         ZDOTC
      // EXTERNAL DZASUM, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, SQRT
      // ..
      // .. Executable Statements ..

      if ( IJOB.NE.2 ) {

         // Apply permutations IPIV to RHS

         CALL ZLASWP( 1, RHS, LDZ, 1, N-1, IPIV, 1 )

         // Solve for L-part choosing RHS either to +1 or -1.

         PMONE = -CONE
         DO 10 J = 1, N - 1
            BP = RHS( J ) + CONE
            BM = RHS( J ) - CONE
            SPLUS = ONE

            // Look-ahead for L- part RHS(1:N-1) = +-1
            // SPLUS and SMIN computed more efficiently than in BSOLVE[1].

            SPLUS = SPLUS + DBLE( ZDOTC( N-J, Z( J+1, J ), 1, Z( J+1, J ), 1 ) )
            SMINU = DBLE( ZDOTC( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 ) )
            SPLUS = SPLUS*DBLE( RHS( J ) )
            if ( SPLUS.GT.SMINU ) {
               RHS( J ) = BP
            } else if ( SMINU.GT.SPLUS ) {
               RHS( J ) = BM
            } else {

               // In this case the updating sums are equal and we can
               // choose RHS(J) +1 or -1. The first time this happens we
               // choose -1, thereafter +1. This is a simple way to get
               // good estimates of matrices like Byers well-known example
               // (see [1]). (Not done in BSOLVE.)

               RHS( J ) = RHS( J ) + PMONE
               PMONE = CONE
            }

            // Compute the remaining r.h.s.

            TEMP = -RHS( J )
            CALL ZAXPY( N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 )
   10    CONTINUE

         // Solve for U- part, lockahead for RHS(N) = +-1. This is not done
         // In BSOLVE and will hopefully give us a better estimate because
         // any ill-conditioning of the original matrix is transferred to U
         // and not to L. U(N, N) is an approximation to sigma_min(LU).

         CALL ZCOPY( N-1, RHS, 1, WORK, 1 )
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
         IF( SPLUS.GT.SMINU ) CALL ZCOPY( N, WORK, 1, RHS, 1 )

         // Apply the permutations JPIV to the computed solution (RHS)

         CALL ZLASWP( 1, RHS, LDZ, 1, N-1, JPIV, -1 )

         // Compute the sum of squares

         CALL ZLASSQ( N, RHS, 1, RDSCAL, RDSUM )
         RETURN
      }

      // ENTRY IJOB = 2

      // Compute approximate nullvector XM of Z

      CALL ZGECON( 'I', N, Z, LDZ, ONE, RTEMP, WORK, RWORK, INFO )
      CALL ZCOPY( N, WORK( N+1 ), 1, XM, 1 )

      // Compute RHS

      CALL ZLASWP( 1, XM, LDZ, 1, N-1, IPIV, -1 )
      TEMP = CONE / SQRT( ZDOTC( N, XM, 1, XM, 1 ) )
      CALL ZSCAL( N, TEMP, XM, 1 )
      CALL ZCOPY( N, XM, 1, XP, 1 )
      CALL ZAXPY( N, CONE, RHS, 1, XP, 1 )
      CALL ZAXPY( N, -CONE, XM, 1, RHS, 1 )
      CALL ZGESC2( N, Z, LDZ, RHS, IPIV, JPIV, SCALE )
      CALL ZGESC2( N, Z, LDZ, XP, IPIV, JPIV, SCALE )
      IF( DZASUM( N, XP, 1 ).GT.DZASUM( N, RHS, 1 ) ) CALL ZCOPY( N, XP, 1, RHS, 1 )

      // Compute the sum of squares

      CALL ZLASSQ( N, RHS, 1, RDSCAL, RDSUM )
      RETURN

      // End of ZLATDF

      }
