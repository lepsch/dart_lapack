      void dchkbb(NSIZES, MVAL, NVAL, NWDTHS, KK, NTYPES, DOTYPE, NRHS, ISEED, THRESH, NOUNIT, A, LDA, AB, LDAB, BD, BE, Q, LDQ, P, LDP, C, LDC, CC, WORK, LWORK, RESULT, INFO ) {

// -- LAPACK test routine (input) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDAB, LDC, LDP, LDQ, LWORK, NOUNIT, NRHS, NSIZES, NTYPES, NWDTHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), MVAL( * ), NVAL( * );
      double             A( LDA, * ), AB( LDAB, * ), BD( * ), BE( * ), C( LDC, * ), CC( LDC, * ), P( LDP, * ), Q( LDQ, * ), RESULT( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                MAXTYP;
      const              MAXTYP = 15 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KL, KMAX, KU, M, MMAX, MNMAX, MNMIN, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDT01, DBDT02, DGBBRD, DLACPY, DLAHD2, DLASET, DLASUM, DLATMR, DLATMS, DORT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 ];
      // ..
      // .. Executable Statements ..

      // Check for errors

      NTESTT = 0;
      INFO = 0;

      // Important constants

      BADMM = false;
      BADNN = false;
      MMAX = 1;
      NMAX = 1;
      MNMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         MMAX = max( MMAX, MVAL( J ) );
         if( MVAL( J ) < 0 ) BADMM = true;
         NMAX = max( NMAX, NVAL( J ) );
         if( NVAL( J ) < 0 ) BADNN = true;
         MNMAX = max( MNMAX, min( MVAL( J ), NVAL( J ) ) );
      } // 10

      BADNNB = false;
      KMAX = 0;
      for (J = 1; J <= NWDTHS; J++) { // 20
         KMAX = max( KMAX, KK( J ) );
         if( KK( J ) < 0 ) BADNNB = true;
      } // 20

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADMM ) {
         INFO = -2;
      } else if ( BADNN ) {
         INFO = -3;
      } else if ( NWDTHS < 0 ) {
         INFO = -4;
      } else if ( BADNNB ) {
         INFO = -5;
      } else if ( NTYPES < 0 ) {
         INFO = -6;
      } else if ( NRHS < 0 ) {
         INFO = -8;
      } else if ( LDA < NMAX ) {
         INFO = -13;
      } else if ( LDAB < 2*KMAX+1 ) {
         INFO = -15;
      } else if ( LDQ < NMAX ) {
         INFO = -19;
      } else if ( LDP < NMAX ) {
         INFO = -21;
      } else if ( LDC < NMAX ) {
         INFO = -23;
      } else if ( ( max( LDA, NMAX )+1 )*NMAX > LWORK ) {
         INFO = -26;
      }

      if ( INFO != 0 ) {
         xerbla('DCHKBB', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      // Loop over sizes, widths, types

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 160
         M = MVAL( JSIZE );
         N = NVAL( JSIZE );
         MNMIN = min( M, N );
         AMNINV = ONE / (max( 1, M, N )).toDouble();

         for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) { // 150
            K = KK( JWIDTH );
            if (K >= M && K >= N) GO TO 150;
            KL = max( 0, min( M-1, K ) );
            KU = max( 0, min( N-1, K ) );

            if ( NSIZES != 1 ) {
               MTYPES = min( MAXTYP, NTYPES );
            } else {
               MTYPES = min( MAXTYP+1, NTYPES );
            }

            for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 140
               if( !DOTYPE( JTYPE ) ) GO TO 140;
               NMATS = NMATS + 1;
               NTEST = 0;

               for (J = 1; J <= 4; J++) { // 30
                  IOLDSD[J] = ISEED( J );
               } // 30

               // Compute "A".

               // Control parameters:

                   // KMAGN  KMODE        KTYPE
               // =1  O(1)   clustered 1  zero
               // =2  large  clustered 2  identity
               // =3  small  exponential  (none)
               // =4         arithmetic   diagonal, (w/ singular values)
               // =5         random log   (none)
               // =6         random       nonhermitian, w/ singular values
               // =7                      (none)
               // =8                      (none)
               // =9                      random nonhermitian

               if (MTYPES > MAXTYP) GO TO 90;

               ITYPE = KTYPE( JTYPE );
               IMODE = KMODE( JTYPE );

               // Compute norm

               GO TO ( 40, 50, 60 )KMAGN( JTYPE );

               } // 40
               ANORM = ONE;
               GO TO 70;

               } // 50
               ANORM = ( RTOVFL*ULP )*AMNINV;
               GO TO 70;

               } // 60
               ANORM = RTUNFL*max( M, N )*ULPINV;
               GO TO 70;

               } // 70

               dlaset('Full', LDA, N, ZERO, ZERO, A, LDA );
               dlaset('Full', LDAB, N, ZERO, ZERO, AB, LDAB );
               IINFO = 0;
               COND = ULPINV;

               // Special Matrices -- Identity & Jordan block

                  // Zero

               if ( ITYPE == 1 ) {
                  IINFO = 0;

               } else if ( ITYPE == 2 ) {

                  // Identity

                  for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                     A[JCOL, JCOL] = ANORM;
                  } // 80

               } else if ( ITYPE == 4 ) {

                  // Diagonal Matrix, singular values specified

                  dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( M+1 ), IINFO );

               } else if ( ITYPE == 6 ) {

                  // Nonhermitian, singular values specified

                  dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, KL, KU, 'N', A, LDA, WORK( M+1 ), IINFO );

               } else if ( ITYPE == 9 ) {

                  // Nonhermitian, random entries

                  dlatmr(M, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, KL, KU, ZERO, ANORM, 'N', A, LDA, IDUMMA, IINFO );

               } else {

                  IINFO = 1;
               }

               // Generate Right-Hand Side

               dlatmr(M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( M+1 ), 1, ONE, WORK( 2*M+1 ), 1, ONE, 'N', IDUMMA, M, NRHS, ZERO, ONE, 'NO', C, LDC, IDUMMA, IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               } // 90

               // Copy A to band storage.

               for (J = 1; J <= N; J++) { // 110
                  for (I = max( 1, J-KU ); I <= min( M, J+KL ); I++) { // 100
                     AB[KU+1+I-J, J] = A( I, J );
                  } // 100
               } // 110

               // Copy C

               dlacpy('Full', M, NRHS, C, LDC, CC, LDC );

               // Call DGBBRD to compute B, Q and P, and to update C.

               dgbbrd('B', M, N, NRHS, KL, KU, AB, LDAB, BD, BE, Q, LDQ, P, LDP, CC, LDC, WORK, IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DGBBRD', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[1] = ULPINV;
                     GO TO 120;
                  }
               }

               // Test 1:  Check the decomposition A := Q * B * P'
                    // 2:  Check the orthogonality of Q
                    // 3:  Check the orthogonality of P
                    // 4:  Check the computation of Q' * C

               dbdt01(M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RESULT( 1 ) );
               dort01('Columns', M, M, Q, LDQ, WORK, LWORK, RESULT( 2 ) );
               dort01('Rows', N, N, P, LDP, WORK, LWORK, RESULT( 3 ) );
               dbdt02(M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RESULT( 4 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 4;
               } // 120
               NTESTT = NTESTT + NTEST;

               // Print out tests which fail.

               for (JR = 1; JR <= NTEST; JR++) { // 130
                  if ( RESULT( JR ) >= THRESH ) {
                     if (NERRS == 0) dlahd2( NOUNIT, 'DBB' );
                     NERRS = NERRS + 1;
                     WRITE( NOUNIT, FMT = 9998 )M, N, K, IOLDSD, JTYPE, JR, RESULT( JR );
                  }
               } // 130

            } // 140
         } // 150
      } // 160

      // Summary

      dlasum('DBB', NOUNIT, NERRS, NTESTT );
      return;

 9999 FORMAT( ' DCHKBB: ', A, ' returned INFO=', I5, '.', / 9X, 'M=', I5, ' N=', I5, ' K=', I5, ', JTYPE=', I5, ', ISEED=(', 3( I5, ',' ), I5, ')' );
 9998 FORMAT( ' M =', I4, ' N=', I4, ', K=', I3, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 );
      }