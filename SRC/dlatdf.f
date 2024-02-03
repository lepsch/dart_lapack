      SUBROUTINE DLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IJOB, LDZ, N;
      double             RDSCAL, RDSUM;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      double             RHS( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                MAXDIM;
      const              MAXDIM = 8 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, K;
      double             BM, BP, PMONE, SMINU, SPLUS, TEMP;
      // ..
      // .. Local Arrays ..
      int                IWORK( MAXDIM );
      double             WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGECON, DGESC2, DLASSQ, DLASWP, DSCAL
      // ..
      // .. External Functions ..
      double             DASUM, DDOT;
      // EXTERNAL DASUM, DDOT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( IJOB.NE.2 ) {

         // Apply permutations IPIV to RHS

         dlaswp(1, RHS, LDZ, 1, N-1, IPIV, 1 );

         // Solve for L-part choosing RHS either to +1 or -1.

         PMONE = -ONE

         for (J = 1; J <= N - 1; J++) { // 10
            BP = RHS( J ) + ONE
            BM = RHS( J ) - ONE
            SPLUS = ONE

            // Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and
            // SMIN computed more efficiently than in BSOLVE [1].

            SPLUS = SPLUS + DDOT( N-J, Z( J+1, J ), 1, Z( J+1, J ), 1 )
            SMINU = DDOT( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 )
            SPLUS = SPLUS*RHS( J )
            if ( SPLUS.GT.SMINU ) {
               RHS( J ) = BP
            } else if ( SMINU.GT.SPLUS ) {
               RHS( J ) = BM
            } else {

               // In this case the updating sums are equal and we can
               // choose RHS(J) +1 or -1. The first time this happens
               // we choose -1, thereafter +1. This is a simple way to
               // get good estimates of matrices like Byers well-known
               // example (see [1]). (Not done in BSOLVE.)

               RHS( J ) = RHS( J ) + PMONE
               PMONE = ONE
            }

            // Compute the remaining r.h.s.

            TEMP = -RHS( J )
            daxpy(N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 );

         } // 10

         // Solve for U-part, look-ahead for RHS(N) = +-1. This is not done
         // in BSOLVE and will hopefully give us a better estimate because
         // any ill-conditioning of the original matrix is transferred to U
         // and not to L. U(N, N) is an approximation to sigma_min(LU).

         dcopy(N-1, RHS, 1, XP, 1 );
         XP( N ) = RHS( N ) + ONE
         RHS( N ) = RHS( N ) - ONE
         SPLUS = ZERO
         SMINU = ZERO
         DO 30 I = N, 1, -1
            TEMP = ONE / Z( I, I )
            XP( I ) = XP( I )*TEMP
            RHS( I ) = RHS( I )*TEMP
            for (K = I + 1; K <= N; K++) { // 20
               XP( I ) = XP( I ) - XP( K )*( Z( I, K )*TEMP )
               RHS( I ) = RHS( I ) - RHS( K )*( Z( I, K )*TEMP )
            } // 20
            SPLUS = SPLUS + ABS( XP( I ) )
            SMINU = SMINU + ABS( RHS( I ) )
         } // 30
         IF( SPLUS.GT.SMINU ) CALL DCOPY( N, XP, 1, RHS, 1 )

         // Apply the permutations JPIV to the computed solution (RHS)

         dlaswp(1, RHS, LDZ, 1, N-1, JPIV, -1 );

         // Compute the sum of squares

         dlassq(N, RHS, 1, RDSCAL, RDSUM );

      } else {

         // IJOB = 2, Compute approximate nullvector XM of Z

         dgecon('I', N, Z, LDZ, ONE, TEMP, WORK, IWORK, INFO );
         dcopy(N, WORK( N+1 ), 1, XM, 1 );

         // Compute RHS

         dlaswp(1, XM, LDZ, 1, N-1, IPIV, -1 );
         TEMP = ONE / SQRT( DDOT( N, XM, 1, XM, 1 ) )
         dscal(N, TEMP, XM, 1 );
         dcopy(N, XM, 1, XP, 1 );
         daxpy(N, ONE, RHS, 1, XP, 1 );
         daxpy(N, -ONE, XM, 1, RHS, 1 );
         dgesc2(N, Z, LDZ, RHS, IPIV, JPIV, TEMP );
         dgesc2(N, Z, LDZ, XP, IPIV, JPIV, TEMP );
         IF( DASUM( N, XP, 1 ).GT.DASUM( N, RHS, 1 ) ) CALL DCOPY( N, XP, 1, RHS, 1 )

         // Compute the sum of squares

         dlassq(N, RHS, 1, RDSCAL, RDSUM );

      }

      RETURN

      // End of DLATDF

      }
