      void claed0(final int QSIZ, final int N, final int D, final int E, final Matrix<double> Q_, final int LDQ, final int QSTORE, final int LDQS, final Array<double> RWORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final Q = Q_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDQ, LDQS, N, QSIZ;
      int                IWORK( * );
      double               D( * ), E( * ), RWORK( * );
      Complex            Q( LDQ, * ), QSTORE( LDQS, * );
      // ..

// =====================================================================

// Warning:      N could be as big as QSIZ!

      // .. Parameters ..
      double               TWO;
      const              TWO = 2.0 ;
      int                CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM, IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM, J, K, LGN, LL, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1, SPM2, SUBMAT, SUBPBS, TLVLS;
      double               TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACRM, CLAED7, SCOPY, SSTEQR, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, REAL

      // Test the input parameters.

      INFO = 0;

      // IF( ICOMPQ < 0 || ICOMPQ > 2 ) THEN
      //    INFO = -1
      // ELSE IF( ( ICOMPQ == 1 ) && ( QSIZ < max( 0, N ) ) )
// $        THEN
      if ( QSIZ < max( 0, N ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDQS < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('CLAED0', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      SMLSIZ = ilaenv( 9, 'CLAED0', ' ', 0, 0, 0, 0 );

      // Determine the size and placement of the submatrices, and save in
      // the leading elements of IWORK.

      IWORK[1] = N;
      SUBPBS = 1;
      TLVLS = 0;
      } // 10
      if ( IWORK( SUBPBS ) > SMLSIZ ) {
         for (J = SUBPBS; J >= 1; J--) { // 20
            IWORK[2*J] = ( IWORK( J )+1 ) / 2;
            IWORK[2*J-1] = IWORK( J ) / 2;
         } // 20
         TLVLS = TLVLS + 1;
         SUBPBS = 2*SUBPBS;
         GO TO 10;
      }
      for (J = 2; J <= SUBPBS; J++) { // 30
         IWORK[J] = IWORK( J ) + IWORK( J-1 );
      } // 30

      // Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
      // using rank-1 modifications (cuts).

      SPM1 = SUBPBS - 1;
      for (I = 1; I <= SPM1; I++) { // 40
         SUBMAT = IWORK( I ) + 1;
         SMM1 = SUBMAT - 1;
         D[SMM1] = D( SMM1 ) - ( E( SMM1 ) ).abs();
         D[SUBMAT] = D( SUBMAT ) - ( E( SMM1 ) ).abs();
      } // 40

      INDXQ = 4*N + 3;

      // Set up workspaces for eigenvalues only/accumulate new vectors
      // routine

      TEMP = LOG( double( N ) ) / LOG( TWO );
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
         IWORK[IPRMPT+I] = 1;
         IWORK[IGIVPT+I] = 1;
      } // 50
      IWORK[IQPTR] = 1;

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
         LL = IQ - 1 + IWORK( IQPTR+CURR );
         ssteqr('I', MATSIZ, D( SUBMAT ), E( SUBMAT ), RWORK( LL ), MATSIZ, RWORK, INFO );
         clacrm(QSIZ, MATSIZ, Q( 1, SUBMAT ), LDQ, RWORK( LL ), MATSIZ, QSTORE( 1, SUBMAT ), LDQS, RWORK( IWREM ) );
         IWORK[IQPTR+CURR+1] = IWORK( IQPTR+CURR ) + MATSIZ**2;
         CURR = CURR + 1;
         if ( INFO > 0 ) {
            INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1;
            return;
         }
         K = 1;
         for (J = SUBMAT; J <= IWORK( I+1 ); J++) { // 60
            IWORK[INDXQ+J] = K;
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
         for (I = 0; I <= SPM2; I += 2) { // 90
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
      // into an eigensystem of size MATSIZ.  CLAED7 handles the case
      // when the eigenvectors of a full or band Hermitian matrix (which
      // was reduced to tridiagonal form) are desired.

      // I am free to use Q as a valuable working space until Loop 150.

            claed7(MATSIZ, MSD2, QSIZ, TLVLS, CURLVL, CURPRB, D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS, E( SUBMAT+MSD2-1 ), IWORK( INDXQ+SUBMAT ), RWORK( IQ ), IWORK( IQPTR ), IWORK( IPRMPT ), IWORK( IPERM ), IWORK( IGIVPT ), IWORK( IGIVCL ), RWORK( IGIVNM ), Q( 1, SUBMAT ), RWORK( IWREM ), IWORK( SUBPBS+1 ), INFO );
            if ( INFO > 0 ) {
               INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1;
               return;
            }
            IWORK[I / 2+1] = IWORK( I+2 );
         } // 90
         SUBPBS = SUBPBS / 2;
         CURLVL = CURLVL + 1;
         GO TO 80;
      }

      // end while

      // Re-merge the eigenvalues/vectors which were deflated at the final
      // merge step.

      for (I = 1; I <= N; I++) { // 100
         J = IWORK( INDXQ+I );
         RWORK[I] = D( J );
         ccopy(QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 );
      } // 100
      scopy(N, RWORK, 1, D, 1 );

      }
