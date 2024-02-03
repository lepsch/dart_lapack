      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      COMPLEX*16         Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, THREE;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) ;
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
      // EXTERNAL DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA, ZLASET, ZLASR, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( LSAME( COMPZ, 'N' ) ) {
         ICOMPZ = 0
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ICOMPZ = 1
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ICOMPZ = 2
      } else {
         ICOMPZ = -1
      }
      if ( ICOMPZ.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZSTEQR', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         IF( ICOMPZ.EQ.2 ) Z( 1, 1 ) = CONE
         RETURN
      }

      // Determine the unit roundoff and over/underflow thresholds.

      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2

      // Compute the eigenvalues and eigenvectors of the tridiagonal
      // matrix.

      IF( ICOMPZ.EQ.2 ) CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )

      NMAXIT = N*MAXIT
      JTOT = 0

      // Determine where the matrix splits and choose QL or QR iteration
      // for each block, according to whether top or bottom diagonal
      // element is smaller.

      L1 = 1
      NM1 = N - 1

      } // 10
      IF( L1.GT.N ) GO TO 160       IF( L1.GT.1 ) E( L1-1 ) = ZERO
      if ( L1.LE.NM1 ) {
         for (M = L1; M <= NM1; M++) { // 20
            TST = ABS( E( M ) )
            if ( TST.EQ.ZERO ) GO TO 30             IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ 1 ) ) ) )*EPS ) {
               E( M ) = ZERO
               GO TO 30
            }
         } // 20
      }
      M = N

      } // 30
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L ) GO TO 10

      // Scale submatrix in rows and columns L to LEND

      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO ) GO TO 10
      if ( ANORM.GT.SSFMAX ) {
         ISCALE = 1
         dlascl('G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )          CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO );
      } else if ( ANORM.LT.SSFMIN ) {
         ISCALE = 2
         dlascl('G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO )          CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO );
      }

      // Choose between QL and QR iteration

      if ( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) {
         LEND = LSV
         L = LENDSV
      }

      if ( LEND.GT.L ) {

         // QL Iteration

         // Look for small subdiagonal element.

         } // 40
         if ( L.NE.LEND ) {
            LENDM1 = LEND - 1
            for (M = L; M <= LENDM1; M++) { // 50
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ SAFMIN )GO TO 60
            } // 50
         }

         M = LEND

         } // 60
         IF( M.LT.LEND ) E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L ) GO TO 80

         // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
         // to compute its eigensystem.

         if ( M.EQ.L+1 ) {
            if ( ICOMPZ.GT.0 ) {
               dlaev2(D( L ), E( L ), D( L+1 ), RT1, RT2, C, S );
               WORK( L ) = C
               WORK( N-1+L ) = S
               zlasr('R', 'V', 'B', N, 2, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ );
            } else {
               dlae2(D( L ), E( L ), D( L+1 ), RT1, RT2 );
            }
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND ) GO TO 40
            GO TO 140
         }

         IF( JTOT.EQ.NMAXIT ) GO TO 140
         JTOT = JTOT + 1

         // Form shift.

         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )

         S = ONE
         C = ONE
         P = ZERO

         // Inner loop

         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            dlartg(G, F, C, S, R );
            IF( I.NE.M-1 ) E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B

            // If eigenvectors are desired, then save rotations.

            if ( ICOMPZ.GT.0 ) {
               WORK( I ) = C
               WORK( N-1+I ) = -S
            }

         } // 70

         // If eigenvectors are desired, then apply saved rotations.

         if ( ICOMPZ.GT.0 ) {
            MM = M - L + 1
            zlasr('R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ );
         }

         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40

         // Eigenvalue found.

         } // 80
         D( L ) = P

         L = L + 1
         IF( L.LE.LEND ) GO TO 40
         GO TO 140

      } else {

         // QR Iteration

         // Look for small superdiagonal element.

         } // 90
         if ( L.NE.LEND ) {
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ SAFMIN )GO TO 110
            } // 100
         }

         M = LEND

         } // 110
         IF( M.GT.LEND ) E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L ) GO TO 130

         // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
         // to compute its eigensystem.

         if ( M.EQ.L-1 ) {
            if ( ICOMPZ.GT.0 ) {
               dlaev2(D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S );
               WORK( M ) = C
               WORK( N-1+M ) = S
               zlasr('R', 'V', 'F', N, 2, WORK( M ), WORK( N-1+M ), Z( 1, L-1 ), LDZ );
            } else {
               dlae2(D( L-1 ), E( L-1 ), D( L ), RT1, RT2 );
            }
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND ) GO TO 90
            GO TO 140
         }

         IF( JTOT.EQ.NMAXIT ) GO TO 140
         JTOT = JTOT + 1

         // Form shift.

         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )

         S = ONE
         C = ONE
         P = ZERO

         // Inner loop

         LM1 = L - 1
         for (I = M; I <= LM1; I++) { // 120
            F = S*E( I )
            B = C*E( I )
            dlartg(G, F, C, S, R );
            IF( I.NE.M ) E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B

            // If eigenvectors are desired, then save rotations.

            if ( ICOMPZ.GT.0 ) {
               WORK( I ) = C
               WORK( N-1+I ) = S
            }

         } // 120

         // If eigenvectors are desired, then apply saved rotations.

         if ( ICOMPZ.GT.0 ) {
            MM = L - M + 1
            zlasr('R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), Z( 1, M ), LDZ );
         }

         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90

         // Eigenvalue found.

         } // 130
         D( L ) = P

         L = L - 1
         IF( L.GE.LEND ) GO TO 90
         GO TO 140

      }

      // Undo scaling if necessary

      } // 140
      if ( ISCALE.EQ.1 ) {
         dlascl('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )          CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO );
      } else if ( ISCALE.EQ.2 ) {
         dlascl('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )          CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO );
      }

      // Check for no convergence to an eigenvalue after a total
      // of N*MAXIT iterations.

      if ( JTOT.EQ.NMAXIT ) {
         for (I = 1; I <= N - 1; I++) { // 150
            IF( E( I ).NE.ZERO ) INFO = INFO + 1
         } // 150
         RETURN
      }
      GO TO 10

      // Order eigenvalues and eigenvectors.

      } // 160
      if ( ICOMPZ.EQ.0 ) {

         // Use Quick Sort

         dlasrt('I', N, D, INFO );

      } else {

         // Use Selection Sort to minimize swaps of eigenvectors

         for (II = 2; II <= N; II++) { // 180
            I = II - 1
            K = I
            P = D( I )
            for (J = II; J <= N; J++) { // 170
               if ( D( J ).LT.P ) {
                  K = J
                  P = D( J )
               }
            } // 170
            if ( K.NE.I ) {
               D( K ) = D( I )
               D( I ) = P
               zswap(N, Z( 1, I ), 1, Z( 1, K ), 1 );
            }
         } // 180
      }
      RETURN

      // End of ZSTEQR

      }
