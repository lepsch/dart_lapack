      void zchkgg(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, TSTDIF, THRSHN, NOUNIT, A, LDA, B, H, T, S1, S2, P1, P2, U, LDU, V, Q, Z, ALPHA1, BETA1, ALPHA3, BETA3, EVECTL, EVECTR, WORK, LWORK, RWORK, LLWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTDIF;
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH, THRSHN;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * ), LLWORK( * );
      int                ISEED( 4 ), NN( * );
      double             RESULT( 15 ), RWORK( * );
      Complex         A( LDA, * ), ALPHA1( * ), ALPHA3( * ), B( LDA, * ), BETA1( * ), BETA3( * ), EVECTL( LDU, * ), EVECTR( LDU, * ), H( LDA, * ), P1( LDA, * ), P2( LDA, * ), Q( LDU, * ), S1( LDA, * ), S2( LDA, * ), T( LDA, * ), U( LDU, * ), V( LDU, * ), WORK( * ), Z( LDU, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 26 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      int                I1, IADD, IINFO, IN, J, JC, JR, JSIZE, JTYPE, LWKOPT, MTYPES, N, N1, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ANORM, BNORM, SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV;
      Complex         CTEMP;
      // ..
      // .. Local Arrays ..
      bool               LASIGN( MAXTYP ), LBSIGN( MAXTYP );
      int                IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ), KATYPE( MAXTYP ), KAZERO( MAXTYP ), KBMAGN( MAXTYP ), KBTYPE( MAXTYP ), KBZERO( MAXTYP ), KCLASS( MAXTYP ), KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 );
      double             DUMMA( 4 ), RMAGN( 0: 3 );
      Complex         CDUMMA( 4 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      //- Complex         ZLARND;
      // EXTERNAL DLAMCH, ZLANGE, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASUM, XERBLA, ZGEQR2, ZGET51, ZGET52, ZGGHRD, ZHGEQZ, ZLACPY, ZLARFG, ZLASET, ZLATM4, ZTGEVC, ZUNM2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN, SIGN
      // ..
      // .. Data statements ..
      const KCLASS = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,];
      const KZ1 = [ 0, 1, 2, 1, 3, 3 ];
      const KZ2 = [ 0, 0, 1, 2, 1, 1 ];
      const KADD = [ 0, 0, 0, 0, 3, 2 ];
      const KATYPE = [ 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0 ];
      const KBTYPE = [ 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0 ];
      const KAZERO = [ 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1 ];
      const KBZERO = [ 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1 ];
      const KAMAGN = [ 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1 ];
      const KBMAGN = [ 1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1 ];
      const KTRIAN = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,];
      const LASIGN = [ false, false, false, false, false, false, true , false , true, true, false, false, true, true, true, false , true , false, false, false, true, true, true, true, true, false ];
      const LBSIGN = [ false, false, false, false, false, false, false, true , false, false, true, true, false, false, true , false , true , false, false, false, false, false, false, false, false, false,];
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0;

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
      } // 10

      LWKOPT = max( 2*NMAX*NMAX, 4*NMAX, 1 );

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADNN ) {
         INFO = -2;
      } else if ( NTYPES < 0 ) {
         INFO = -3;
      } else if ( THRESH < ZERO ) {
         INFO = -6;
      } else if ( LDA <= 1 || LDA < NMAX ) {
         INFO = -10;
      } else if ( LDU <= 1 || LDU < NMAX ) {
         INFO = -19;
      } else if ( LWKOPT > LWORK ) {
         INFO = -30;
      }

      if ( INFO != 0 ) {
         xerbla('ZCHKGG', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      SAFMIN = dlamch( 'Safe minimum' );
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );
      SAFMIN = SAFMIN / ULP;
      SAFMAX = ONE / SAFMIN;
      ULPINV = ONE / ULP;

      // The values RMAGN(2:3) depend on N, see below.

      RMAGN[0] = ZERO;
      RMAGN[1] = ONE;

      // Loop over sizes, types

      NTESTT = 0;
      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 240
         N = NN( JSIZE );
         N1 = max( 1, N );
         RMAGN[2] = SAFMAX*ULP / N1.toDouble();
         RMAGN[3] = SAFMIN*ULPINV*N1;

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 230
            if( !DOTYPE( JTYPE ) ) GO TO 230;
            NMATS = NMATS + 1;
            NTEST = 0;

            // Save ISEED in case of an error.

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            // Initialize RESULT

            for (J = 1; J <= 15; J++) { // 30
               RESULT[J] = ZERO;
            } // 30

            // Compute A and B

            // Description of control parameters:

            // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
                    // =3 means random.
            // KATYPE: the "type" to be passed to ZLATM4 for computing A.
            // KAZERO: the pattern of zeros on the diagonal for A:
                    // =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
                    // =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
                    // =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
                    // non-zero entries.)
            // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
                    // =2: large, =3: small.
            // LASIGN: true if the diagonal elements of A are to be
                    // multiplied by a random magnitude 1 number.
            // KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN:  used to implement KAMAGN and KBMAGN.

            if (MTYPES > MAXTYP) GO TO 110;
            IINFO = 0;
            if ( KCLASS( JTYPE ) < 3 ) {

               // Generate A (w/o rotation)

               if ( ( KATYPE( JTYPE ) ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) zlaset( 'Full', N, N, CZERO, CZERO, A, LDA );
               } else {
                  IN = N;
               }
               zlatm4(KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ), KZ2( KAZERO( JTYPE ) ), LASIGN( JTYPE ), RMAGN( KAMAGN( JTYPE ) ), ULP, RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 4, ISEED, A, LDA );
               IADD = KADD( KAZERO( JTYPE ) );
               if (IADD > 0 && IADD <= N) A( IADD, IADD ) = RMAGN( KAMAGN( JTYPE ) );

               // Generate B (w/o rotation)

               if ( ( KBTYPE( JTYPE ) ).abs() == 3 ) {
                  IN = 2*( ( N-1 ) / 2 ) + 1;
                  if (IN != N) zlaset( 'Full', N, N, CZERO, CZERO, B, LDA );
               } else {
                  IN = N;
               }
               zlatm4(KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ), KZ2( KBZERO( JTYPE ) ), LBSIGN( JTYPE ), RMAGN( KBMAGN( JTYPE ) ), ONE, RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 4, ISEED, B, LDA );
               IADD = KADD( KBZERO( JTYPE ) );
               if (IADD != 0) B( IADD, IADD ) = RMAGN( KBMAGN( JTYPE ) );

               if ( KCLASS( JTYPE ) == 2 && N > 0 ) {

                  // Include rotations

                  // Generate U, V as Householder transformations times a
                  // diagonal matrix.  (Note that ZLARFG makes U(j,j) and
                  // V(j,j) real.)

                  for (JC = 1; JC <= N - 1; JC++) { // 50
                     for (JR = JC; JR <= N; JR++) { // 40
                        U[JR][JC] = ZLARND( 3, ISEED );
                        V[JR][JC] = ZLARND( 3, ISEED );
                     } // 40
                     zlarfg(N+1-JC, U( JC, JC ), U( JC+1, JC ), 1, WORK( JC ) );
                     WORK[2*N+JC] = sign( ONE, (U( JC, JC )).toDouble() );
                     U[JC][JC] = CONE;
                     zlarfg(N+1-JC, V( JC, JC ), V( JC+1, JC ), 1, WORK( N+JC ) );
                     WORK[3*N+JC] = sign( ONE, (V( JC, JC )).toDouble() );
                     V[JC][JC] = CONE;
                  } // 50
                  CTEMP = ZLARND( 3, ISEED );
                  U[N][N] = CONE;
                  WORK[N] = CZERO;
                  WORK[3*N] = CTEMP / ( CTEMP ).abs();
                  CTEMP = ZLARND( 3, ISEED );
                  V[N][N] = CONE;
                  WORK[2*N] = CZERO;
                  WORK[4*N] = CTEMP / ( CTEMP ).abs();

                  // Apply the diagonal matrices

                  for (JC = 1; JC <= N; JC++) { // 70
                     for (JR = 1; JR <= N; JR++) { // 60
                        A[JR][JC] = WORK( 2*N+JR )* DCONJG( WORK( 3*N+JC ) )* A( JR, JC )                         B( JR, JC ) = WORK( 2*N+JR )* DCONJG( WORK( 3*N+JC ) )* B( JR, JC );
                     } // 60
                  } // 70
                  CALL ZUNM2R( 'L', 'N', N, N, N-1, U, LDU, WORK, A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL ZUNM2R( 'R', 'C', N, N, N-1, V, LDU, WORK( N+1 ), A, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL ZUNM2R( 'L', 'N', N, N, N-1, U, LDU, WORK, B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
                  CALL ZUNM2R( 'R', 'C', N, N, N-1, V, LDU, WORK( N+1 ), B, LDA, WORK( 2*N+1 ), IINFO )                   IF( IINFO != 0 ) GO TO 100;
               }
            } else {

               // Random matrices

               for (JC = 1; JC <= N; JC++) { // 90
                  for (JR = 1; JR <= N; JR++) { // 80
                     A[JR][JC] = RMAGN( KAMAGN( JTYPE ) )* ZLARND( 4, ISEED )                      B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )* ZLARND( 4, ISEED );
                  } // 80
               } // 90
            }

            ANORM = ZLANGE( '1', N, N, A, LDA, RWORK );
            BNORM = ZLANGE( '1', N, N, B, LDA, RWORK );

            } // 100

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 110

            // Call ZGEQR2, ZUNM2R, and ZGGHRD to compute H, T, U, and V

            zlacpy(' ', N, N, A, LDA, H, LDA );
            zlacpy(' ', N, N, B, LDA, T, LDA );
            NTEST = 1;
            RESULT[1] = ULPINV;

            zgeqr2(N, N, T, LDA, WORK, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZGEQR2', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zunm2r('L', 'C', N, N, N, T, LDA, WORK, H, LDA, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZUNM2R', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zlaset('Full', N, N, CZERO, CONE, U, LDU );
            zunm2r('R', 'N', N, N, N, T, LDA, WORK, U, LDU, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZUNM2R', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zgghrd('V', 'I', N, 1, N, H, LDA, T, LDA, U, LDU, V, LDU, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZGGHRD', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }
            NTEST = 4;

            // Do tests 1--4

            zget51(1, N, A, LDA, H, LDA, U, LDU, V, LDU, WORK, RWORK, RESULT( 1 ) );
            zget51(1, N, B, LDA, T, LDA, U, LDU, V, LDU, WORK, RWORK, RESULT( 2 ) );
            zget51(3, N, B, LDA, T, LDA, U, LDU, U, LDU, WORK, RWORK, RESULT( 3 ) );
            zget51(3, N, B, LDA, T, LDA, V, LDU, V, LDU, WORK, RWORK, RESULT( 4 ) );

            // Call ZHGEQZ to compute S1, P1, S2, P2, Q, and Z, do tests.

            // Compute T1 and UZ

            // Eigenvalues only

            zlacpy(' ', N, N, H, LDA, S2, LDA );
            zlacpy(' ', N, N, T, LDA, P2, LDA );
            NTEST = 5;
            RESULT[5] = ULPINV;

            zhgeqz('E', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHA3, BETA3, Q, LDU, Z, LDU, WORK, LWORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHGEQZ(E)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            // Eigenvalues and Full Schur Form

            zlacpy(' ', N, N, H, LDA, S2, LDA );
            zlacpy(' ', N, N, T, LDA, P2, LDA );

            zhgeqz('S', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHA1, BETA1, Q, LDU, Z, LDU, WORK, LWORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHGEQZ(S)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            // Eigenvalues, Schur Form, and Schur Vectors

            zlacpy(' ', N, N, H, LDA, S1, LDA );
            zlacpy(' ', N, N, T, LDA, P1, LDA );

            zhgeqz('S', 'I', 'I', N, 1, N, S1, LDA, P1, LDA, ALPHA1, BETA1, Q, LDU, Z, LDU, WORK, LWORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZHGEQZ(V)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            NTEST = 8;

            // Do Tests 5--8

            zget51(1, N, H, LDA, S1, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 5 ) );
            zget51(1, N, T, LDA, P1, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 6 ) );
            zget51(3, N, T, LDA, P1, LDA, Q, LDU, Q, LDU, WORK, RWORK, RESULT( 7 ) );
            zget51(3, N, T, LDA, P1, LDA, Z, LDU, Z, LDU, WORK, RWORK, RESULT( 8 ) );

            // Compute the Left and Right Eigenvectors of (S1,P1)

            // 9: Compute the left eigenvector Matrix without
               // back transforming:

            NTEST = 9;
            RESULT[9] = ULPINV;

            // To test "SELECT" option, compute half of the eigenvectors
            // in one call, and half in another

            I1 = N / 2;
            for (J = 1; J <= I1; J++) { // 120
               LLWORK[J] = true;
            } // 120
            for (J = I1 + 1; J <= N; J++) { // 130
               LLWORK[J] = false;
            } // 130

            ztgevc('L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTGEVC(L,S1)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            I1 = IN;
            for (J = 1; J <= I1; J++) { // 140
               LLWORK[J] = false;
            } // 140
            for (J = I1 + 1; J <= N; J++) { // 150
               LLWORK[J] = true;
            } // 150

            ztgevc('L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL( 1, I1+1 ), LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTGEVC(L,S2)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zget52( true , N, S1, LDA, P1, LDA, EVECTL, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) );
            RESULT[9] = DUMMA( 1 );
            if ( DUMMA( 2 ) > THRSHN ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZTGEVC(HOWMNY=S)', DUMMA( 2 ), N, JTYPE, IOLDSD;
            }

            // 10: Compute the left eigenvector Matrix with
                // back transforming:

            NTEST = 10;
            RESULT[10] = ULPINV;
            zlacpy('F', N, N, Q, LDU, EVECTL, LDU );
            ztgevc('L', 'B', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU, CDUMMA, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTGEVC(L,B)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zget52( true , N, H, LDA, T, LDA, EVECTL, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) );
            RESULT[10] = DUMMA( 1 );
            if ( DUMMA( 2 ) > THRSHN ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'ZTGEVC(HOWMNY=B)', DUMMA( 2 ), N, JTYPE, IOLDSD;
            }

            // 11: Compute the right eigenvector Matrix without
                // back transforming:

            NTEST = 11;
            RESULT[11] = ULPINV;

            // To test "SELECT" option, compute half of the eigenvectors
            // in one call, and half in another

            I1 = N / 2;
            for (J = 1; J <= I1; J++) { // 160
               LLWORK[J] = true;
            } // 160
            for (J = I1 + 1; J <= N; J++) { // 170
               LLWORK[J] = false;
            } // 170

            ztgevc('R', 'S', LLWORK, N, S1, LDA, P1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTGEVC(R,S1)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            I1 = IN;
            for (J = 1; J <= I1; J++) { // 180
               LLWORK[J] = false;
            } // 180
            for (J = I1 + 1; J <= N; J++) { // 190
               LLWORK[J] = true;
            } // 190

            ztgevc('R', 'S', LLWORK, N, S1, LDA, P1, LDA, CDUMMA, LDU, EVECTR( 1, I1+1 ), LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTGEVC(R,S2)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zget52( false , N, S1, LDA, P1, LDA, EVECTR, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) );
            RESULT[11] = DUMMA( 1 );
            if ( DUMMA( 2 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZTGEVC(HOWMNY=S)', DUMMA( 2 ), N, JTYPE, IOLDSD;
            }

            // 12: Compute the right eigenvector Matrix with
                // back transforming:

            NTEST = 12;
            RESULT[12] = ULPINV;
            zlacpy('F', N, N, Z, LDU, EVECTR, LDU );
            ztgevc('R', 'B', LLWORK, N, S1, LDA, P1, LDA, CDUMMA, LDU, EVECTR, LDU, N, IN, WORK, RWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'ZTGEVC(R,B)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               GO TO 210;
            }

            zget52( false , N, H, LDA, T, LDA, EVECTR, LDU, ALPHA1, BETA1, WORK, RWORK, DUMMA( 1 ) );
            RESULT[12] = DUMMA( 1 );
            if ( DUMMA( 2 ) > THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'ZTGEVC(HOWMNY=B)', DUMMA( 2 ), N, JTYPE, IOLDSD;
            }

            // Tests 13--15 are done only on request

            if ( TSTDIF ) {

               // Do Tests 13--14

               zget51(2, N, S1, LDA, S2, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 13 ) );
               zget51(2, N, P1, LDA, P2, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT( 14 ) );

               // Do Test 15

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 200
                  TEMP1 = max( TEMP1, ABS( ALPHA1( J )-ALPHA3( J ) ) );
                  TEMP2 = max( TEMP2, ABS( BETA1( J )-BETA3( J ) ) );
               } // 200

               TEMP1 = TEMP1 / max( SAFMIN, ULP*max( TEMP1, ANORM ) );
               TEMP2 = TEMP2 / max( SAFMIN, ULP*max( TEMP2, BNORM ) );
               RESULT[15] = max( TEMP1, TEMP2 );
               NTEST = 15;
            } else {
               RESULT[13] = ZERO;
               RESULT[14] = ZERO;
               RESULT[15] = ZERO;
               NTEST = 12;
            }

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 210

            NTESTT = NTESTT + NTEST;

            // Print out tests which fail.

            for (JR = 1; JR <= NTEST; JR++) { // 220
               if ( RESULT( JR ) >= THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9997 )'ZGG';

                     // Matrix types

                     WRITE( NOUNIT, FMT = 9996 );
                     WRITE( NOUNIT, FMT = 9995 );
                     WRITE( NOUNIT, FMT = 9994 )'Unitary';

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9993 )'unitary', '*', 'conjugate transpose', ( '*', J = 1, 10 );

                  }
                  NERRS = NERRS + 1;
                  if ( RESULT( JR ) < 10000.0 ) {
                     WRITE( NOUNIT, FMT = 9992 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  } else {
                     WRITE( NOUNIT, FMT = 9991 )N, JTYPE, IOLDSD, JR, RESULT( JR );
                  }
               }
            } // 220

         } // 230
      } // 240

      // Summary

      dlasum('ZGG', NOUNIT, NERRS, NTESTT );
      return;

 9999 FORMAT( ' ZCHKGG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

 9998 FORMAT( ' ZCHKGG: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

 9997 FORMAT( 1X, A3, ' -- Complex Generalized eigenvalue problem' );

 9996 FORMAT( ' Matrix types (see ZCHKGG for details): ' );

 9995 FORMAT( ' Special Matrices:', 23X, '(J''=transposed Jordan block)', / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ', '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ', 'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I', ')  11=(large*I, small*D)  13=(large*D, large*I)', / '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ', ' 14=(small*D, small*I)', / '  15=(D, reversed D)' );
 9994 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:', / '  16=Transposed Jordan Blocks             19=geometric ', 'alpha, beta=0,1', / '  17=arithm. alpha&beta             ', '      20=arithmetic alpha, beta=0,1', / '  18=clustered ', 'alpha, beta=0,1            21=random alpha, beta=0,1', / ' Large & Small Matrices:', / '  22=(large, small)   ', '23=(small,large)    24=(small,small)    25=(large,large)', / '  26=random O(1) matrices.' );

 9993 FORMAT( / ' Tests performed:   (H is Hessenberg, S is Schur, B, ', 'T, P are triangular,', / 20X, 'U, V, Q, and Z are ', A, ', l and r are the', / 20X, 'appropriate left and right eigenvectors, resp., a is', / 20X, 'alpha, b is beta, and ', A, ' means ', A, '.)', / ' 1 = | A - U H V', A, ' | / ( |A| n ulp )      2 = | B - U T V', A, ' | / ( |B| n ulp )', / ' 3 = | I - UU', A, ' | / ( n ulp )             4 = | I - VV', A, ' | / ( n ulp )', / ' 5 = | H - Q S Z', A, ' | / ( |H| n ulp )', 6X, '6 = | T - Q P Z', A, ' | / ( |T| n ulp )', / ' 7 = | I - QQ', A, ' | / ( n ulp )             8 = | I - ZZ', A, ' | / ( n ulp )', / ' 9 = max | ( b S - a P )', A, ' l | / const.  10 = max | ( b H - a T )', A, ' l | / const.', / ' 11= max | ( b S - a P ) r | / const.   12 = max | ( b H', ' - a T ) r | / const.', / 1X );

 9992 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 0P, F8.2 );
 9991 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I2, ' is', 1P, D10.3 );
      }
