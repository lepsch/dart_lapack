      SUBROUTINE SLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, WORK, IWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDQ, LDQS, N, QSIZ;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM, IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM, J, K, LGN, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1, SPM2, SUBMAT, SUBPBS, TLVLS;
      REAL               TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SLACPY, SLAED1, SLAED7, SSTEQR, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( ICOMPQ < 0 || ICOMPQ > 2 ) {
         INFO = -1;
      } else if ( ( ICOMPQ == 1 ) && ( QSIZ < max( 0, N ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDQS < max( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SLAED0', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      SMLSIZ = ILAENV( 9, 'SLAED0', ' ', 0, 0, 0, 0 );

      // Determine the size and placement of the submatrices, and save in
      // the leading elements of IWORK.

      IWORK( 1 ) = N;
      SUBPBS = 1;
      TLVLS = 0;
      } // 10
      if ( IWORK( SUBPBS ) > SMLSIZ ) {
         DO 20 J = SUBPBS, 1, -1;
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2;
            IWORK( 2*J-1 ) = IWORK( J ) / 2;
         } // 20
         TLVLS = TLVLS + 1;
         SUBPBS = 2*SUBPBS;
         GO TO 10;
      }
      for (J = 2; J <= SUBPBS; J++) { // 30
         IWORK( J ) = IWORK( J ) + IWORK( J-1 );
      } // 30

      // Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
      // using rank-1 modifications (cuts).

      SPM1 = SUBPBS - 1;
      for (I = 1; I <= SPM1; I++) { // 40
         SUBMAT = IWORK( I ) + 1;
         SMM1 = SUBMAT - 1;
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) );
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) );
      } // 40

      INDXQ = 4*N + 3;
      if ( ICOMPQ != 2 ) {

         // Set up workspaces for eigenvalues only/accumulate new vectors
         // routine

         TEMP = LOG( REAL( N ) ) / LOG( TWO );
         LGN = INT( TEMP );
         if (2**LGN < N) LGN = LGN + 1;
         IF( 2**LGN < N ) LGN = LGN + 1;
         IPRMPT = INDXQ + N + 1;
         IPERM = IPRMPT + N*LGN;
         IQPTR = IPERM + N*LGN;
         IGIVPT = IQPTR + N + 2;
         IGIVCL = IGIVPT + N*LGN;

         IGIVNM = 1;
         IQ = IGIVNM + 2*N*LGN;
         IWREM = IQ + N**2 + 1;

         // Initialize pointers

         for (I = 0; I <= SUBPBS; I++) { // 50
            IWORK( IPRMPT+I ) = 1;
            IWORK( IGIVPT+I ) = 1;
         } // 50
         IWORK( IQPTR ) = 1;
      }

      // Solve each submatrix eigenproblem at the bottom of the divide and
      // conquer tree.

      CURR = 0;
      for (I = 0; I <= SPM1; I++) { // 70
         if ( I == 0 ) {
            SUBMAT = 1;
            MATSIZ = IWORK( 1 );
         } else {
            SUBMAT = IWORK( I ) + 1;
            MATSIZ = IWORK( I+1 ) - IWORK( I );
         }
         if ( ICOMPQ == 2 ) {
            CALL SSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ), Q( SUBMAT, SUBMAT ), LDQ, WORK, INFO )             IF( INFO != 0 ) GO TO 130;
         } else {
            ssteqr('I', MATSIZ, D( SUBMAT ), E( SUBMAT ), WORK( IQ-1+IWORK( IQPTR+CURR ) ), MATSIZ, WORK, INFO );
            if (INFO != 0) GO TO 130;
            if ( ICOMPQ == 1 ) {
               sgemm('N', 'N', QSIZ, MATSIZ, MATSIZ, ONE, Q( 1, SUBMAT ), LDQ, WORK( IQ-1+IWORK( IQPTR+ CURR ) ), MATSIZ, ZERO, QSTORE( 1, SUBMAT ), LDQS );
            }
            IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2;
            CURR = CURR + 1;
         }
         K = 1;
         for (J = SUBMAT; J <= IWORK( I+1 ); J++) { // 60
            IWORK( INDXQ+J ) = K;
            K = K + 1;
         } // 60
      } // 70

      // Successively merge eigensystems of adjacent submatrices
      // into eigensystem for the corresponding larger matrix.

      // while ( SUBPBS > 1 )

      CURLVL = 1;
      } // 80
      if ( SUBPBS > 1 ) {
         SPM2 = SUBPBS - 2;
         DO 90 I = 0, SPM2, 2;
            if ( I == 0 ) {
               SUBMAT = 1;
               MATSIZ = IWORK( 2 );
               MSD2 = IWORK( 1 );
               CURPRB = 0;
            } else {
               SUBMAT = IWORK( I ) + 1;
               MATSIZ = IWORK( I+2 ) - IWORK( I );
               MSD2 = MATSIZ / 2;
               CURPRB = CURPRB + 1;
            }

      // Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
      // into an eigensystem of size MATSIZ.
      // SLAED1 is used only for the full eigensystem of a tridiagonal
      // matrix.
      // SLAED7 handles the cases in which eigenvalues only or eigenvalues
      // and eigenvectors of a full symmetric matrix (which was reduced to
      // tridiagonal form) are desired.

            if ( ICOMPQ == 2 ) {
               slaed1(MATSIZ, D( SUBMAT ), Q( SUBMAT, SUBMAT ), LDQ, IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ), MSD2, WORK, IWORK( SUBPBS+1 ), INFO );
            } else {
               slaed7(ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB, D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS, IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ), MSD2, WORK( IQ ), IWORK( IQPTR ), IWORK( IPRMPT ), IWORK( IPERM ), IWORK( IGIVPT ), IWORK( IGIVCL ), WORK( IGIVNM ), WORK( IWREM ), IWORK( SUBPBS+1 ), INFO );
            }
            if (INFO != 0) GO TO 130;
            IWORK( I / 2+1 ) = IWORK( I+2 );
         } // 90
         SUBPBS = SUBPBS / 2;
         CURLVL = CURLVL + 1;
         GO TO 80;
      }

      // end while

      // Re-merge the eigenvalues/vectors which were deflated at the final
      // merge step.

      if ( ICOMPQ == 1 ) {
         for (I = 1; I <= N; I++) { // 100
            J = IWORK( INDXQ+I );
            WORK( I ) = D( J );
            scopy(QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 );
         } // 100
         scopy(N, WORK, 1, D, 1 );
      } else if ( ICOMPQ == 2 ) {
         for (I = 1; I <= N; I++) { // 110
            J = IWORK( INDXQ+I );
            WORK( I ) = D( J );
            scopy(N, Q( 1, J ), 1, WORK( N*I+1 ), 1 );
         } // 110
         scopy(N, WORK, 1, D, 1 );
         slacpy('A', N, N, WORK( N+1 ), N, Q, LDQ );
      } else {
         for (I = 1; I <= N; I++) { // 120
            J = IWORK( INDXQ+I );
            WORK( I ) = D( J );
         } // 120
         scopy(N, WORK, 1, D, 1 );
      }
      GO TO 140;

      } // 130
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1;

      } // 140
      return;

      // End of SLAED0

      }
