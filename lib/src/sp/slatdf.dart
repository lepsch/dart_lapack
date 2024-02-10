      void slatdf(IJOB, N, final Matrix<double> Z, final int LDZ, RHS, RDSUM, RDSCAL, final Array<int> IPIV, final int JPIV) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IJOB, LDZ, N;
      double               RDSCAL, RDSUM;
      int                IPIV( * ), JPIV( * );
      double               RHS( * ), Z( LDZ, * );
      // ..

      int                MAXDIM;
      const              MAXDIM = 8 ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J, K;
      double               BM, BP, PMONE, SMINU, SPLUS, TEMP;
      int                IWORK( MAXDIM );
      double               WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGECON, SGESC2, SLASSQ, SLASWP, SSCAL
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SDOT;
      // EXTERNAL SASUM, SDOT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT

      if ( IJOB != 2 ) {

         // Apply permutations IPIV to RHS

         slaswp(1, RHS, LDZ, 1, N-1, IPIV, 1 );

         // Solve for L-part choosing RHS either to +1 or -1.

         PMONE = -ONE;

         for (J = 1; J <= N - 1; J++) { // 10
            BP = RHS( J ) + ONE;
            BM = RHS( J ) - ONE;
            SPLUS = ONE;

            // Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and
            // SMIN computed more efficiently than in BSOLVE [1].

            SPLUS = SPLUS + SDOT( N-J, Z( J+1, J ), 1, Z( J+1, J ), 1 );
            SMINU = SDOT( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 );
            SPLUS = SPLUS*RHS( J );
            if ( SPLUS > SMINU ) {
               RHS[J] = BP;
            } else if ( SMINU > SPLUS ) {
               RHS[J] = BM;
            } else {

               // In this case the updating sums are equal and we can
               // choose RHS(J) +1 or -1. The first time this happens
               // we choose -1, thereafter +1. This is a simple way to
               // get good estimates of matrices like Byers well-known
               // example (see [1]). (Not done in BSOLVE.)

               RHS[J] = RHS( J ) + PMONE;
               PMONE = ONE;
            }

            // Compute the remaining r.h.s.

            TEMP = -RHS( J );
            saxpy(N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 );

         } // 10

         // Solve for U-part, look-ahead for RHS(N) = +-1. This is not done
         // in BSOLVE and will hopefully give us a better estimate because
         // any ill-conditioning of the original matrix is transferred to U
         // and not to L. U(N, N) is an approximation to sigma_min(LU).

         scopy(N-1, RHS, 1, XP, 1 );
         XP[N] = RHS( N ) + ONE;
         RHS[N] = RHS( N ) - ONE;
         SPLUS = ZERO;
         SMINU = ZERO;
         for (I = N; I >= 1; I--) { // 30
            TEMP = ONE / Z( I, I );
            XP[I] = XP( I )*TEMP;
            RHS[I] = RHS( I )*TEMP;
            for (K = I + 1; K <= N; K++) { // 20
               XP[I] = XP( I ) - XP( K )*( Z( I, K )*TEMP );
               RHS[I] = RHS( I ) - RHS( K )*( Z( I, K )*TEMP );
            } // 20
            SPLUS = SPLUS + ( XP( I ) ).abs();
            SMINU = SMINU + ( RHS( I ) ).abs();
         } // 30
         if (SPLUS > SMINU) scopy( N, XP, 1, RHS, 1 );

         // Apply the permutations JPIV to the computed solution (RHS)

         slaswp(1, RHS, LDZ, 1, N-1, JPIV, -1 );

         // Compute the sum of squares

         slassq(N, RHS, 1, RDSCAL, RDSUM );

      } else {

         // IJOB = 2, Compute approximate nullvector XM of Z

         sgecon('I', N, Z, LDZ, ONE, TEMP, WORK, IWORK, INFO );
         scopy(N, WORK( N+1 ), 1, XM, 1 );

         // Compute RHS

         slaswp(1, XM, LDZ, 1, N-1, IPIV, -1 );
         TEMP = ONE / sqrt( SDOT( N, XM, 1, XM, 1 ) );
         sscal(N, TEMP, XM, 1 );
         scopy(N, XM, 1, XP, 1 );
         saxpy(N, ONE, RHS, 1, XP, 1 );
         saxpy(N, -ONE, XM, 1, RHS, 1 );
         sgesc2(N, Z, LDZ, RHS, IPIV, JPIV, TEMP );
         sgesc2(N, Z, LDZ, XP, IPIV, JPIV, TEMP );
         if( SASUM( N, XP, 1 ) > SASUM( N, RHS, 1 ) ) scopy( N, XP, 1, RHS, 1 );

         // Compute the sum of squares

         slassq(N, RHS, 1, RDSCAL, RDSUM );

      }

      }
