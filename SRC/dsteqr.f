      void dsteqr(COMPZ, N, D, E, Z, LDZ, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, THREE;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0 ;
      int                MAXIT;
      const              MAXIT = 30 ;
      // ..
      // .. Local Scalars ..
      int                I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, NM1, NMAXIT;
      double             ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANST, DLAPY2;
      // EXTERNAL LSAME, DLAMCH, DLANST, DLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR, DLASRT, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( LSAME( COMPZ, 'N' ) ) {
         ICOMPZ = 0;
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ICOMPZ = 1;
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ICOMPZ = 2;
      } else {
         ICOMPZ = -1;
      }
      if ( ICOMPZ < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ( LDZ < 1 ) || ( ICOMPZ > 0 && LDZ < max( 1, N ) ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DSTEQR', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if (ICOMPZ == 2) Z( 1, 1 ) = ONE;
         return;
      }

      // Determine the unit roundoff and over/underflow thresholds.

      EPS = DLAMCH( 'E' );
      EPS2 = EPS**2;
      SAFMIN = DLAMCH( 'S' );
      SAFMAX = ONE / SAFMIN;
      SSFMAX = sqrt( SAFMAX ) / THREE;
      SSFMIN = sqrt( SAFMIN ) / EPS2;

      // Compute the eigenvalues and eigenvectors of the tridiagonal
      // matrix.

      if (ICOMPZ == 2) CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ );

      NMAXIT = N*MAXIT;
      JTOT = 0;

      // Determine where the matrix splits and choose QL or QR iteration
      // for each block, according to whether top or bottom diagonal
      // element is smaller.

      L1 = 1;
      NM1 = N - 1;

      } // 10
      if (L1 > N) GO TO 160;
      IF( L1 > 1 ) E( L1-1 ) = ZERO;
      if ( L1 <= NM1 ) {
         for (M = L1; M <= NM1; M++) { // 20
            TST = ABS( E( M ) );
            if ( TST == ZERO ) GO TO 30;
            IF( TST <= ( sqrt( ABS( D( M ) ) )*sqrt( ABS( D( M+ 1 ) ) ) )*EPS ) {
               E( M ) = ZERO;
               GO TO 30;
            }
         } // 20
      }
      M = N;

      } // 30
      L = L1;
      LSV = L;
      LEND = M;
      LENDSV = LEND;
      L1 = M + 1;
      if (LEND == L) GO TO 10;

      // Scale submatrix in rows and columns L to LEND

      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) );
      ISCALE = 0;
      if (ANORM == ZERO) GO TO 10;
      if ( ANORM > SSFMAX ) {
         ISCALE = 1;
         dlascl('G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO );
         dlascl('G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO );
      } else if ( ANORM < SSFMIN ) {
         ISCALE = 2;
         dlascl('G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO );
         dlascl('G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO );
      }

      // Choose between QL and QR iteration

      if ( ABS( D( LEND ) ) < ABS( D( L ) ) ) {
         LEND = LSV;
         L = LENDSV;
      }

      if ( LEND > L ) {

         // QL Iteration

         // Look for small subdiagonal element.

         } // 40
         if ( L != LEND ) {
            LENDM1 = LEND - 1;
            for (M = L; M <= LENDM1; M++) { // 50
               TST = ABS( E( M ) )**2;
               if( TST <= ( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ SAFMIN )GO TO 60;
            } // 50
         }

         M = LEND;

         } // 60
         if (M < LEND) E( M ) = ZERO;
         P = D( L );
         if (M == L) GO TO 80;

         // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
         // to compute its eigensystem.

         if ( M == L+1 ) {
            if ( ICOMPZ > 0 ) {
               dlaev2(D( L ), E( L ), D( L+1 ), RT1, RT2, C, S );
               WORK( L ) = C;
               WORK( N-1+L ) = S;
               dlasr('R', 'V', 'B', N, 2, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ );
            } else {
               dlae2(D( L ), E( L ), D( L+1 ), RT1, RT2 );
            }
            D( L ) = RT1;
            D( L+1 ) = RT2;
            E( L ) = ZERO;
            L = L + 2;
            if (L <= LEND) GO TO 40;
            GO TO 140;
         }

         if (JTOT == NMAXIT) GO TO 140;
         JTOT = JTOT + 1;

         // Form shift.

         G = ( D( L+1 )-P ) / ( TWO*E( L ) );
         R = DLAPY2( G, ONE );
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) );

         S = ONE;
         C = ONE;
         P = ZERO;

         // Inner loop

         MM1 = M - 1;
         DO 70 I = MM1, L, -1;
            F = S*E( I );
            B = C*E( I );
            dlartg(G, F, C, S, R );
            if (I != M-1) E( I+1 ) = R;
            G = D( I+1 ) - P;
            R = ( D( I )-G )*S + TWO*C*B;
            P = S*R;
            D( I+1 ) = G + P;
            G = C*R - B;

            // If eigenvectors are desired, then save rotations.

            if ( ICOMPZ > 0 ) {
               WORK( I ) = C;
               WORK( N-1+I ) = -S;
            }

         } // 70

         // If eigenvectors are desired, then apply saved rotations.

         if ( ICOMPZ > 0 ) {
            MM = M - L + 1;
            dlasr('R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ );
         }

         D( L ) = D( L ) - P;
         E( L ) = G;
         GO TO 40;

         // Eigenvalue found.

         } // 80
         D( L ) = P;

         L = L + 1;
         if (L <= LEND) GO TO 40;
         GO TO 140;

      } else {

         // QR Iteration

         // Look for small superdiagonal element.

         } // 90
         if ( L != LEND ) {
            LENDP1 = LEND + 1;
            DO 100 M = L, LENDP1, -1;
               TST = ABS( E( M-1 ) )**2;
               if( TST <= ( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ SAFMIN )GO TO 110;
            } // 100
         }

         M = LEND;

         } // 110
         if (M > LEND) E( M-1 ) = ZERO;
         P = D( L );
         if (M == L) GO TO 130;

         // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
         // to compute its eigensystem.

         if ( M == L-1 ) {
            if ( ICOMPZ > 0 ) {
               dlaev2(D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S );
               WORK( M ) = C;
               WORK( N-1+M ) = S;
               dlasr('R', 'V', 'F', N, 2, WORK( M ), WORK( N-1+M ), Z( 1, L-1 ), LDZ );
            } else {
               dlae2(D( L-1 ), E( L-1 ), D( L ), RT1, RT2 );
            }
            D( L-1 ) = RT1;
            D( L ) = RT2;
            E( L-1 ) = ZERO;
            L = L - 2;
            if (L >= LEND) GO TO 90;
            GO TO 140;
         }

         if (JTOT == NMAXIT) GO TO 140;
         JTOT = JTOT + 1;

         // Form shift.

         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) );
         R = DLAPY2( G, ONE );
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) );

         S = ONE;
         C = ONE;
         P = ZERO;

         // Inner loop

         LM1 = L - 1;
         for (I = M; I <= LM1; I++) { // 120
            F = S*E( I );
            B = C*E( I );
            dlartg(G, F, C, S, R );
            if (I != M) E( I-1 ) = R;
            G = D( I ) - P;
            R = ( D( I+1 )-G )*S + TWO*C*B;
            P = S*R;
            D( I ) = G + P;
            G = C*R - B;

            // If eigenvectors are desired, then save rotations.

            if ( ICOMPZ > 0 ) {
               WORK( I ) = C;
               WORK( N-1+I ) = S;
            }

         } // 120

         // If eigenvectors are desired, then apply saved rotations.

         if ( ICOMPZ > 0 ) {
            MM = L - M + 1;
            dlasr('R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), Z( 1, M ), LDZ );
         }

         D( L ) = D( L ) - P;
         E( LM1 ) = G;
         GO TO 90;

         // Eigenvalue found.

         } // 130
         D( L ) = P;

         L = L - 1;
         if (L >= LEND) GO TO 90;
         GO TO 140;

      }

      // Undo scaling if necessary

      } // 140
      if ( ISCALE == 1 ) {
         dlascl('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO );
         dlascl('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO );
      } else if ( ISCALE == 2 ) {
         dlascl('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO );
         dlascl('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO );
      }

      // Check for no convergence to an eigenvalue after a total
      // of N*MAXIT iterations.

      if (JTOT < NMAXIT) GO TO 10;
      for (I = 1; I <= N - 1; I++) { // 150
         if( E( I ) != ZERO ) INFO = INFO + 1;
      } // 150
      GO TO 190;

      // Order eigenvalues and eigenvectors.

      } // 160
      if ( ICOMPZ == 0 ) {

         // Use Quick Sort

         dlasrt('I', N, D, INFO );

      } else {

         // Use Selection Sort to minimize swaps of eigenvectors

         for (II = 2; II <= N; II++) { // 180
            I = II - 1;
            K = I;
            P = D( I );
            for (J = II; J <= N; J++) { // 170
               if ( D( J ) < P ) {
                  K = J;
                  P = D( J );
               }
            } // 170
            if ( K != I ) {
               D( K ) = D( I );
               D( I ) = P;
               dswap(N, Z( 1, I ), 1, Z( 1, K ), 1 );
            }
         } // 180
      }

      } // 190
      return;

      // End of DSTEQR

      }
