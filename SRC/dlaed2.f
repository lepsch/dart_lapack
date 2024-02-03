      SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMBDA, W, Q2, INDX, INDXC, INDXP, COLTYP, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDQ, N, N1;
      double             RHO;
      // ..
      // .. Array Arguments ..
      int                COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ), INDXQ( * );
      double             D( * ), DLAMBDA( * ), Q( LDQ, * ), Q2( * ), W( * ), Z( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             MONE, ZERO, ONE, TWO, EIGHT;
      const              MONE = -1.0, ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0 ;
      // ..
      // .. Local Arrays ..
      int                CTOT( 4 ), PSM( 4 );
      // ..
      // .. Local Scalars ..
      int                CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1, N2, NJ, PJ;
      double             C, EPS, S, T, TAU, TOL;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DLAPY2;
      // EXTERNAL IDAMAX, DLAMCH, DLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( N < 0 ) {
         INFO = -2;
      } else if ( LDQ < MAX( 1, N ) ) {
         INFO = -6;
      } else if ( MIN( 1, ( N / 2 ) ) > N1 || ( N / 2 ) < N1 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('DLAED2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      N2 = N - N1;
      N1P1 = N1 + 1;

      if ( RHO < ZERO ) {
         dscal(N2, MONE, Z( N1P1 ), 1 );
      }

      // Normalize z so that norm(z) = 1.  Since z is the concatenation of
      // two normalized vectors, norm2(z) = sqrt(2).

      T = ONE / SQRT( TWO );
      dscal(N, T, Z, 1 );

      // RHO = ABS( norm(z)**2 * RHO )

      RHO = ABS( TWO*RHO );

      // Sort the eigenvalues into increasing order

      for (I = N1P1; I <= N; I++) { // 10
         INDXQ( I ) = INDXQ( I ) + N1;
      } // 10

      // re-integrate the deflated parts from the last pass

      for (I = 1; I <= N; I++) { // 20
         DLAMBDA( I ) = D( INDXQ( I ) );
      } // 20
      dlamrg(N1, N2, DLAMBDA, 1, 1, INDXC );
      for (I = 1; I <= N; I++) { // 30
         INDX( I ) = INDXQ( INDXC( I ) );
      } // 30

      // Calculate the allowable deflation tolerance

      IMAX = IDAMAX( N, Z, 1 );
      JMAX = IDAMAX( N, D, 1 );
      EPS = DLAMCH( 'Epsilon' );
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) );

      // If the rank-1 modifier is small enough, no more needs to be done
      // except to reorganize Q so that its columns correspond with the
      // elements in D.

      if ( RHO*ABS( Z( IMAX ) ) <= TOL ) {
         K = 0;
         IQ2 = 1;
         for (J = 1; J <= N; J++) { // 40
            I = INDX( J );
            dcopy(N, Q( 1, I ), 1, Q2( IQ2 ), 1 );
            DLAMBDA( J ) = D( I );
            IQ2 = IQ2 + N;
         } // 40
         dlacpy('A', N, N, Q2, N, Q, LDQ );
         dcopy(N, DLAMBDA, 1, D, 1 );
         GO TO 190;
      }

      // If there are multiple eigenvalues then the problem deflates.  Here
      // the number of equal eigenvalues are found.  As each equal
      // eigenvalue is found, an elementary reflector is computed to rotate
      // the corresponding eigensubspace so that the corresponding
      // components of Z are zero in this new basis.

      for (I = 1; I <= N1; I++) { // 50
         COLTYP( I ) = 1;
      } // 50
      for (I = N1P1; I <= N; I++) { // 60
         COLTYP( I ) = 3;
      } // 60


      K = 0;
      K2 = N + 1;
      for (J = 1; J <= N; J++) { // 70
         NJ = INDX( J );
         if ( RHO*ABS( Z( NJ ) ) <= TOL ) {

            // Deflate due to small z component.

            K2 = K2 - 1;
            COLTYP( NJ ) = 4;
            INDXP( K2 ) = NJ;
            if (J == N) GO TO 100;
         } else {
            PJ = NJ;
            GO TO 80;
         }
      } // 70
      } // 80
      J = J + 1;
      NJ = INDX( J );
      if (J > N) GO TO 100;
      if ( RHO*ABS( Z( NJ ) ) <= TOL ) {

         // Deflate due to small z component.

         K2 = K2 - 1;
         COLTYP( NJ ) = 4;
         INDXP( K2 ) = NJ;
      } else {

         // Check if eigenvalues are close enough to allow deflation.

         S = Z( PJ );
         C = Z( NJ );

         // Find sqrt(a**2+b**2) without overflow or
         // destructive underflow.

         TAU = DLAPY2( C, S );
         T = D( NJ ) - D( PJ );
         C = C / TAU;
         S = -S / TAU;
         if ( ABS( T*C*S ) <= TOL ) {

            // Deflation is possible.

            Z( NJ ) = TAU;
            Z( PJ ) = ZERO;
            if( COLTYP( NJ ) != COLTYP( PJ ) ) COLTYP( NJ ) = 2;
            COLTYP( PJ ) = 4;
            drot(N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S );
            T = D( PJ )*C**2 + D( NJ )*S**2;
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2;
            D( PJ ) = T;
            K2 = K2 - 1;
            I = 1;
            } // 90
            if ( K2+I <= N ) {
               if ( D( PJ ) < D( INDXP( K2+I ) ) ) {
                  INDXP( K2+I-1 ) = INDXP( K2+I );
                  INDXP( K2+I ) = PJ;
                  I = I + 1;
                  GO TO 90;
               } else {
                  INDXP( K2+I-1 ) = PJ;
               }
            } else {
               INDXP( K2+I-1 ) = PJ;
            }
            PJ = NJ;
         } else {
            K = K + 1;
            DLAMBDA( K ) = D( PJ );
            W( K ) = Z( PJ );
            INDXP( K ) = PJ;
            PJ = NJ;
         }
      }
      GO TO 80;
      } // 100

      // Record the last eigenvalue.

      K = K + 1;
      DLAMBDA( K ) = D( PJ );
      W( K ) = Z( PJ );
      INDXP( K ) = PJ;

      // Count up the total number of the various types of columns, then
      // form a permutation which positions the four column types into
      // four uniform groups (although one or more of these groups may be
      // empty).

      for (J = 1; J <= 4; J++) { // 110
         CTOT( J ) = 0;
      } // 110
      for (J = 1; J <= N; J++) { // 120
         CT = COLTYP( J );
         CTOT( CT ) = CTOT( CT ) + 1;
      } // 120

      // PSM(*) = Position in SubMatrix (of types 1 through 4)

      PSM( 1 ) = 1;
      PSM( 2 ) = 1 + CTOT( 1 );
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 );
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 );
      K = N - CTOT( 4 );

      // Fill out the INDXC array so that the permutation which it induces
      // will place all type-1 columns first, all type-2 columns next,
      // then all type-3's, and finally all type-4's.

      for (J = 1; J <= N; J++) { // 130
         JS = INDXP( J );
         CT = COLTYP( JS );
         INDX( PSM( CT ) ) = JS;
         INDXC( PSM( CT ) ) = J;
         PSM( CT ) = PSM( CT ) + 1;
      } // 130

      // Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
      // and Q2 respectively.  The eigenvalues/vectors which were not
      // deflated go into the first K slots of DLAMBDA and Q2 respectively,
      // while those which were deflated go into the last N - K slots.

      I = 1;
      IQ1 = 1;
      IQ2 = 1 + ( CTOT( 1 )+CTOT( 2 ) )*N1;
      for (J = 1; J <= CTOT( 1 ); J++) { // 140
         JS = INDX( I );
         dcopy(N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 );
         Z( I ) = D( JS );
         I = I + 1;
         IQ1 = IQ1 + N1;
      } // 140

      for (J = 1; J <= CTOT( 2 ); J++) { // 150
         JS = INDX( I );
         dcopy(N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 );
         dcopy(N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 );
         Z( I ) = D( JS );
         I = I + 1;
         IQ1 = IQ1 + N1;
         IQ2 = IQ2 + N2;
      } // 150

      for (J = 1; J <= CTOT( 3 ); J++) { // 160
         JS = INDX( I );
         dcopy(N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 );
         Z( I ) = D( JS );
         I = I + 1;
         IQ2 = IQ2 + N2;
      } // 160

      IQ1 = IQ2;
      for (J = 1; J <= CTOT( 4 ); J++) { // 170
         JS = INDX( I );
         dcopy(N, Q( 1, JS ), 1, Q2( IQ2 ), 1 );
         IQ2 = IQ2 + N;
         Z( I ) = D( JS );
         I = I + 1;
      } // 170

      // The deflated eigenvalues and their corresponding vectors go back
      // into the last N - K slots of D and Q respectively.

      if ( K < N ) {
         dlacpy('A', N, CTOT( 4 ), Q2( IQ1 ), N, Q( 1, K+1 ), LDQ );
         dcopy(N-K, Z( K+1 ), 1, D( K+1 ), 1 );
      }

      // Copy CTOT into COLTYP for referencing in DLAED3.

      for (J = 1; J <= 4; J++) { // 180
         COLTYP( J ) = CTOT( J );
      } // 180

      } // 190
      return;

      // End of DLAED2

      }
